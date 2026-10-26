# SAM Epigenetics Pipeline: Raw IDAT Analysis of GSE112893

This repository contains the reproducible R scripts used to analyze raw DNA methylation data (Illumina 450k IDAT files) from the GSE112893 cohort, investigating epigenetic dysregulation in Edematous Severe Acute Malnutrition (ESAM).

## 📂 Repository Contents
- `IDAT_Gold_Standard.R`: Main pipeline for reading raw IDATs, Noob background correction, quality control, and `limma` differential methylation analysis.
- `Recover_and_Annotate.R`: Script for extracting significant DMPs, mapping probes to genomic features using Bioconductor annotation packages, and generating gene lists for pathway enrichment.

## 🛠️ Requirements
- R (>= 4.5.0)
- Bioconductor packages: `minfi`, `limma`, `IlluminaHumanMethylation450kanno.ilmn12.hg19`, `IlluminaHumanMethylation450kmanifest`
- CRAN packages: `tidyverse`, `ggplot2`

## 🚀 How to Run
1. Download the raw IDAT files for GSE112893 from the [NCBI GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112893).
2. Place the `.idat` files in your working directory.
3. Open RStudio, set your working directory, and source the scripts:
   ```R
   source("IDAT_Gold_Standard.R")
   source("Recover_and_Annotate.R")
