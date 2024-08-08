
# Failure to rescue orthoebolavirus-exposed macaques with antibody 1C3 coincides with rapid emergence of resistance mutations

> In brief, five macaques were challenged with SUDV Gulu ([P2 from clinical isolate](https://doi.org/10.1172/jci.insight.159090)). Three macaques were treated with antibody 1C3, which [targets the "chalice" of the EBOV glycoprotein](https://doi.org/10.1016/j.cell.2022.02.023); the two remaining macaques were untreated controls. We sequenced rRNA-depleted RNA from four tissues from all five macaques: axillary lymph node (AxLN), liver, plasma, and spleen. The resulting sequences were examined for potential resistance mutations to 1C3 treatment.

## RNA and sequencing library prep

Total RNA was extracted from either plasma or homogenized tissue. RNA was cleaned, and DNA was digested, with the Zymo Clean and Concentrator-5 kit. Sequencing libraries were constructed following rRNA depletion and sequenced on an Element Biosciences Aviti sequencer to a target depth of 50 million reads per sample.

## Data processing pipeline

The software and pipeline used are outlined below. For full code and details, please see [`step1-process.sh`](https://github.com/geisbert-lab/cyno-sudv-prot1902012-1/blob/9c219e9fb1589dbfe591ad3cc22216fc7df69421/step1-process.sh).

### Pipeline overview

1. Host-aligning reads were removed by aligning to the *M. fascicularis* genome ([v6.0](https://useast.ensembl.org/Macaca_fascicularis/Info/Index)) using [`bowtie2`](https://doi.org/10.1038/nmeth.1923) with the `--very-sensitive` flag. 
2. Reads that did not align to the cyno genome were aligned to the SUDV Gulu RefSeq ([NC_006432.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_006432.1)) using [`bowtie2`](https://doi.org/10.1038/nmeth.1923).
3. The viral alignment BAM file was sorted via [`samtools sort`](https://www.htslib.org/doc/samtools-sort.html).
4. Indel mutations were assessed with [`lofreq indelqual`](https://doi.org/10.1093/nar/gks918), and the resulting BAM was indexed using [`samtools index`](https://www.htslib.org/doc/samtools-index.html).
5. SUDV Gulu genome coverage was calculated using the BAM file and [`samtools depth`](https://www.htslib.org/doc/samtools-depth.html).
6. SNVs were called and quantified with [`lofreq call`](https://doi.org/10.1093/nar/gks918). Only SNVs with >100 aligned reads were reported.

### Software used

| Package    | Version |
| :--------  | :------ |
| `bowtie2`  | 2.4.2   |
| `samtools` | 1.18    |
| `lofreq`   | 2.1.3.1 |

## SNV analysis

SNV analyses were performed in R as outlined below. For full code and details, please see `step2-analyze.r`.

### Analysis overview

1. SNVs were filtered (frequency >0.001) and annotated with [`Biostrings`](https://doi.org/doi:10.18129/B9.bioc.Biostrings) to classify noncoding, synonymous, and nonsynonymous mutations. Special handling was included for the [SUDV Gulu GP transcriptional editing site](https://www.ncbi.nlm.nih.gov/nuccore/NC_006432.1#feature_NC_006432.1_misc_feature_0). 
2. [SUDV Gulu GP coverage](https://github.com/geisbert-lab/cyno-sudv-prot1902012-1/blob/9c219e9fb1589dbfe591ad3cc22216fc7df69421/analysis/coverage-gp.png) was checked to confirm that all nucleotides in all samples had adequate coverage for SNV calling (>100 aligned reads).
3. [GP's SNV profile](https://github.com/geisbert-lab/cyno-sudv-prot1902012-1/blob/9c219e9fb1589dbfe591ad3cc22216fc7df69421/analysis/snvs-gp.png) was plotted for all samples to visualize the overall intrahost variation.
4. Potential resistance mutations were identified by filtering for SNVs that met the following criteria:
	- Nonsynonymous
	- In GP
	- Found at >0.1 frequency in at least two treated NHPs
	- *Not* found at >0.1 frequency in any control sample
5. One SNV, [GP P124L](https://github.com/geisbert-lab/cyno-sudv-prot1902012-1/blob/9c219e9fb1589dbfe591ad3cc22216fc7df69421/analysis/snvs-gpp124L-column.png), met these criteria. The difference in GP P124L frequences between control and treated animal tissues was statistically significant ([Wilcoxon nonparamtreic test p<0.001](https://github.com/geisbert-lab/cyno-sudv-prot1902012-1/blob/9c219e9fb1589dbfe591ad3cc22216fc7df69421/analysis/snvs-gpp124L-violin.png)).

### Software used

| Package     | Version |
| :---------- | :------ |
| `R`         | 4.0.2   |
| `tidyverse` | 2.0.0   |
| `Biostrings`| 2.56.0  |
| `ggpubr`    | 0.6.0   |

