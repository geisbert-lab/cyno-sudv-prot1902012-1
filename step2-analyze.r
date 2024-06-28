#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(argparse))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

# helper variables
anot <- read.csv("pipeline/genomes/sudv-gulu-NC006432.1.csv")
cols.cond <- c(Control="darkgrey", Treated="#525252")

# helper functions
read.vcf <- function(id) {
  paste0("data/", id, "/snvs.vcf") %>%
    read.delim(comment.char="#", header=FALSE,
               col.names=c("Segment", "NT.position", 
                           "NT.ID", "NT.ref", "NT.alt", 
                           "QUAL", "FILTER", "INFO")) %>%
    # quality filter
    filter(FILTER=="PASS") %>%
    # extract frequency, format NT ID, and add sample ID
    mutate(Frequency=str_extract(INFO, "(?<=AF=)[0-9\\.]+"),
           Frequency=as.numeric(Frequency),
           NT.ID=paste0(NT.position, "-", NT.ref, "-", NT.alt),
           ID=id) %>%
    # format table and return
    select(ID, NT.position, NT.ID, Frequency, NT.ref, NT.alt)
}
read.cov <- function(id) {
  paste0("data/", id, "/coverage.tsv") %>%
    read.delim(col.names=c("Segment", "NT.position", "Depth")) %>%
    # add sample ID
    mutate(ID=id) %>%
    # format table and return
    select(ID, NT.position, Depth)
}
annotate.snvs <- function(nt.id, annotation=anot, stutter=6883) {
  # extract nucleotide info
  nt.pos <- as.integer(str_extract(nt.id, ("^[0-9]+")))
  nt.ref <- str_extract(nt.id, "(?<=-)[ACGT]+(?=-)")
  nt.alt <- str_extract(nt.id, "[ACGT]+$")
  
  # init return data frame
  res.df <- data.frame(NT.position=nt.pos,
                       NT.ID=nt.id,
                       NT.ref=nt.ref,
                       NT.alt=nt.alt,
                       Type=NA,
                       Gene=NA,
                       AA.ID=NA,
                       AA.position=NA,
                       AA.ref=NA,
                       AA.alt=NA)
  
  # is the mutation coding or noncoding?
  # if it's noncoding, return
  orf <- filter(annotation, Start <= nt.pos, End >= nt.pos)
  if(dim(orf)[1]==0) {
    res.df$Type <- "Noncoding"
    res.df$AA.ID <- paste("Noncoding", nt.id)
    return(res.df)
  }
  
  # edge case check: GP stutter handling
  # if the NT.position is after the stutter site, we need to shift the ORF
  # starting position -1 to keep the # nucleotides and coordinates consistent
  if(orf$Gene=="GP" & nt.pos > stutter) {
    orf$Start <- orf$Start-1
  }
  
  # format reference AA sequence
  aa.ref <- orf$AASeq %>% 
    strsplit("") %>%
    unlist()
  
  # format mutated NT sequence to have NT positions as indices
  aa.alt <- orf$NucSeq %>%
    strsplit("") %>%
    unlist()
  names(aa.alt) <- orf$Start:orf$End
  
  # edge case check: deletions (insertions don't affect reference)
  if(nchar(nt.ref) > 1) {
    # how many nucleotide to delete?
    ndel <- nchar(nt.ref) - 1
    # what are the coordinated of these deleted nucleotides?
    ndel <- (nt.pos+1):(nt.pos+ndel)
    # substitute NTs for NAs
    aa.alt[as.character(ndel)] <- NA
    # re-make NT sequence without NAs
    aa.alt <- na.omit(aa.alt)
  } else { # substitute to mutate NT sequence
    aa.alt[as.character(nt.pos)] <- nt.alt
  }
  
  # translate mutated NT sequence to AA
  aa.alt <- aa.alt %>%
    paste(collapse="") %>%
    DNAString() %>%
    translate() %>% 
    as.character() %>%
    strsplit("") %>%
    unlist()
  
  # check for nonsynonymous changes
  aa.pos <- which(aa.ref != aa.alt)
  
  # if no AA mutations, format return data frame for synonymous
  if(length(aa.pos)==0) {
    res.df$Type <- "Synonymous"
    res.df$Gene <- orf$Gene
    res.df$AA.ID <- paste(orf$Gene, "synonymous", nt.id)
    return(res.df)
  }
  
  # if we're still standing, it's a nonsynonymous mutation
  # quick check to grab only the first AA if a frameshift
  aa.pos <- aa.pos[1]
  aa.ref <- aa.ref[aa.pos]
  aa.alt <- aa.alt[aa.pos]
  # format data frame and return
  res.df$Type <- "Nonsynonymous"
  res.df$Gene <- orf$Gene
  res.df$AA.ID <- paste0(orf$Gene, " ", aa.ref, aa.pos, aa.alt)
  res.df$AA.position <- aa.pos
  res.df$AA.ref <- aa.ref
  res.df$AA.alt <- aa.alt
  return(res.df)
}

## inputs ----------------------------------------------------------------------
# metadata
meta <- read.csv("analysis/samplesheet.csv") %>%
        select(-starts_with("RNA")) %>%
        mutate(NHP=paste("NHP", NHP),
               NHP=factor(NHP, levels=c("NHP 131338", "NHP 1401953",
                                        "NHP 1210733", "NHP 1210855",
                                        "NHP 1306692")),
               Condition=factor(Condition, levels=c("Control", "Treated")))

# read in VCF files and filter/format
snvs <- lapply(meta$ID, read.vcf) %>%
        do.call(rbind, .) %>%
        # remove SNVs with frequency < Illumina error (0.001)
        filter(Frequency > 0.001)

# annotate SNVs, format, and save full dataset
snvs <- snvs$NT.ID %>%
        unique() %>%
        lapply(annotate.snvs) %>%
        do.call(rbind, .) %>% 
        right_join(snvs, 
                   by=c("NT.position", "NT.ID", "NT.ref", "NT.alt")) %>%
        select(ID, NT.position, Frequency, AA.ID, NT.ID, 
               Type, Gene, NT.ref, NT.alt, AA.position, AA.ref, AA.alt) %>%
        dplyr::arrange(NT.position)
write.csv(snvs, "analysis/snvs.csv", na="", row.names=FALSE)

# read in coverage and build profile
prof <- lapply(meta$ID, read.cov) %>%
        do.call(rbind, .)
# collapse and format SNVs
x <- snvs %>%
     group_by(ID, NT.position) %>%
     summarise(Frequency=sum(Frequency),
               .groups="drop")
# merge coverage and SNV frequency
prof <- left_join(prof, x, by=c("ID", "NT.position")) %>%
        # add in zeros
        replace_na(list(Frequency=0)) %>%
        # add in metadata
        left_join(meta, by="ID")
rm(x)

## QC coverage and frequency (whole genome) ------------------------------------
# coverage
prof %>%
  ggplot(aes(NT.position, Depth, col=Condition)) +
  geom_line() +
  scale_color_manual(values=cols.cond) +
  geom_hline(yintercept=100, linetype=2, col="red") +
  scale_y_continuous(limits=c(1, 1.5e6), trans="log10") +
  facet_grid(NHP ~ Tissue) +
  labs(x="Nucleotide position",
       y="Aligned read depth")
ggsave("analysis/coverage-all.png", units="cm", width=25, height=15)

# snvs
prof %>%
  ggplot(aes(NT.position, Frequency, col=Condition)) +
  geom_line() +
  scale_color_manual(values=cols.cond) +
  geom_hline(yintercept=0.5, linetype=2, col="red") +
  facet_grid(NHP ~ Tissue) +
  labs(x="Nucleotide position",
       y="SNV frequency")
ggsave("analysis/snvs-all.png", units="cm", width=25, height=15)

## QC coverage and frequency (GP only) -----------------------------------------
# subset everything to GP
snvs <- filter(snvs, Gene=="GP")
prof <- filter(prof,
               NT.position >= anot[anot$Gene=="GP", "Start"],
               NT.position <= anot[anot$Gene=="GP", "End"])

# coverage
prof %>%
  ggplot(aes(NT.position, Depth, col=Condition)) +
  geom_line() +
  scale_color_manual(values=cols.cond) +
  geom_hline(yintercept=100, linetype=2, col="red") +
  scale_y_continuous(limits=c(1, 1.5e6), trans="log10") +
  facet_grid(NHP ~ Tissue) +
  labs(x="Nucleotide position",
       y="Aligned read depth")
ggsave("analysis/coverage-gp.png", units="cm", width=25, height=15)

# snvs
prof %>%
  ggplot(aes(NT.position, Frequency, col=Condition)) +
  geom_line() +
  scale_color_manual(values=cols.cond) +
  geom_hline(yintercept=0.5, linetype=2, col="red") +
  ylim(0, 1) +
  facet_grid(NHP ~ Tissue) +
  labs(x="Nucleotide position",
       y="SNV frequency")
ggsave("analysis/snvs-gp.png", units="cm", width=25, height=15)

## investigating potentially adaptive SNVs -------------------------------------
# only care about nonsynonymous SNVs > 10% in GP
# also get rid of 7U/8U
target <- filter(snvs, 
                 Frequency > 0.1, 
                 Type=="Nonsynonymous",
                 NT.ID != "6876-T-TA") %>%
          left_join(meta, by="ID")

# filter to get SNVs only found in treated
baseline <- target %>%
            filter(Condition=="Control") %>%
            select(AA.ID) %>%
            unlist()
adaptive <- target %>%
            filter(Condition=="Treated",
                   !(AA.ID %in% baseline)) %>%
            select(ID, AA.ID, Frequency)
# expand to fill in zeros. use SNVs data frame to fill in frequencies <0.1 if present
x <- expand.grid(ID=meta[meta$Condition=="Treated", "ID"],
                 AA.ID=unique(adaptive$AA.ID))
adaptive <- left_join(x, snvs, by=colnames(x)) %>%
            left_join(meta, by="ID") %>%
            replace_na(list(Frequency=0))

# plot adaptive mutation frequencies
adaptive %>%
  mutate(NHP=str_remove(NHP, "^NHP")) %>%
  ggplot(aes(NHP, Frequency, fill=Tissue)) +
  geom_col(position="dodge", col="black") +
  scale_fill_brewer(palette="Paired") +
  geom_hline(yintercept=0.1, col="red", linetype=2) +
  facet_wrap(~AA.ID) +
  ylim(0, 1) +
  labs(x="Treated NHPs",
       y="SNV frequency") +
  theme(legend.position=c(0.15, 0.75),
        legend.background=element_rect(color="black"))
ggsave("analysis/snvs-treated.png",
       units="cm", width=15, height=10)

# GP P124L found in multiple treated NHPs at >10%
adaptive %>%
  filter(AA.ID=="GP P124L") %>%
  mutate(NHP=str_remove(NHP, "^NHP")) %>%
  ggplot(aes(NHP, Frequency, fill=Tissue)) +
  geom_col(position="dodge", col="black") +
  scale_fill_brewer(palette="Paired") +
  geom_hline(yintercept=0.1, col="red", linetype=2) +
  ylim(0, 1) +
  labs(x="Treated NHPs",
       y="GP P124L frequency") +
  theme(legend.position=c(0.15, 0.8),
        legend.background=element_rect(color="black"))
ggsave("analysis/snvs-gpp124L-column.png",
       units="cm", width=10, height=10)

# are the GP P124L frequencies statistically higher in treated than controls?
snvs %>%
  filter(AA.ID=="GP P124L") %>%
  select(ID, Frequency) %>%
  right_join(meta, by="ID") %>%
  # fill in missing zeros
  replace_na(list(Frequency=0)) %>%
  ggplot(aes(Condition, Frequency, fill=Condition)) +
  geom_boxplot(alpha=0.5) +
  geom_jitter(pch=21, height=0, width=0.1) +
  scale_fill_manual(values=cols.cond) +
  ylim(0, 1) +
  ggpubr::stat_compare_means(comparisons=list(c("Control", "Treated")), 
                             label.y=0.6, label="p.format") +
  labs(x=element_blank(),
       y="GP P124L frequency") +
  theme(legend.position="none")
ggsave("analysis/snvs-gpp124L-violin.png",
       units="cm", width=5, height=10)

# what's the coverage at this position? NT 6368
prof %>%
  filter(NT.position==6368) %>%
  mutate(NHP=str_remove(NHP, "^NHP")) %>%
  ggplot(aes(NHP, Depth, group=Tissue, fill=Condition)) +
  geom_col(col="black", position="dodge") +
  scale_fill_manual(values=cols.cond) +
  geom_hline(yintercept=100, col="red", linetype=2) +
  scale_y_continuous(limits=c(1, 1.5e6), trans="log10") +
  labs(x=element_blank(),
       y="GP P124L aligned read depth")
ggsave("analysis/coverage-gpp124L.png",
       units="cm", width=15, height=10)

# write CSV with coverage
adaptive %>%
  filter(AA.ID=="GP P124L") %>%
  mutate(NHP=str_remove(NHP, "^NHP")) %>%
  select(ID, NHP, Tissue, Frequency, NT.position) %>%
  left_join(select(prof, ID, NT.position, Depth),
            by=c("ID", "NT.position")) %>%
  select(-NT.position) %>%
  write.csv("analysis/snvs-gpp124L.csv",
          row.names=FALSE)