# Identifying potential resistence mutatants in SUDV-challenged cynomolgus macaques treated with 1C3

> In brief, five macaques were challenged with SUDV Gulu ([P2 from clinical isolate](https://doi.org/10.1172/jci.insight.159090)). Three macaques were treated with antibody 1C3, which [targets the "chalice" of the EBOV glycoprotein](https://doi.org/10.1016/j.cell.2022.02.023); the two remaining macaques were untreated controls. We sequenced rRNA-depleted RNA from four tissues from all five macaques: axillary lymph node (AxLN), liver, plasma, and spleen. The resulting sequences were examined for potential resistance mutations to 1C3 treatment.

## RNA and sequencing library prep

Total RNA was extracted from either plasma or homogenized tissue. All RNA samples were cleaned, and DNA was digested, with the Zymo Clean and Concentrator-5 kit. Sequencing libraries were constructed following rRNA depletion and sequenced on an Element Biosciences Aviti sequencer to a target depth of 50 million reads per sample.

## Data processing pipeline

Software and pipeline are outlined below. For full code and details, please see `step1-process.sh`

### Pipeline overview

1. Host-aligning reads were removed by aligning to the *M. fascicularis* genome ([v6.0](https://useast.ensembl.org/Macaca_fascicularis/Info/Index)) using `bowtie2` with the `--very-sensitive` flag. 
2. Reads that did not align to the cyno genome were aligned to the SUDV Gulu RefSeq ([NC_006432.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_006432.1)) using `bowtie2`.
3. The viral alignment BAM file was sorted via `samtools sort`.
4. Indel mutations were assessed with `lofreq indelqual`, and the resulting BAM was indexed using `samtools index`.
5. SUDV Gulu genome coverage was calculated using the BAM file and `samtools depth`.
6. SNVs were called and quantified with `lofreq call`.

### Software used

| Package    | Version |
| :--------: | :-----: |
| bowtie2    | 2.4.2   |
| samtools   | 1.18    |
| lofreq     | 2.1.3.1 |

## SNV analysis

...to be continued...

