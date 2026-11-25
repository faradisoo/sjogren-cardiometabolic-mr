# B-cell transcriptomics and Mendelian randomisation in Sjögren's syndrome

Analysis code for:

> Muhammad Agil Wijaya Faradis. Sjögren's syndrome and cardiometabolic risk: a Mendelian randomization study integrating salivary gland transcriptomics. *medRxiv* (preprint), 2026. DOI: [to be added after posting]

The study combines three in silico analyses:

1. **Transcriptomics.** Differential expression in parotid gland tissue (GEO GSE40611: 17 primary Sjögren's syndrome [pSS], 14 non-pSS sicca, 18 controls; Affymetrix HG-U133 Plus 2.0), followed by GO and KEGG enrichment.
2. **MR analysis A.** Genetic liability to Sjögren's syndrome → seven cardiometabolic traits (LDL, HDL, triglycerides, CRP, systolic and diastolic blood pressure, coronary heart disease).
3. **MR analysis B.** Blood cis-eQTLs of B-cell genes up-regulated in pSS (eQTLGen) → Sjögren's syndrome and the same cardiometabolic traits, with colocalisation for findings that pass Bonferroni correction.

## Repository contents

| File | Purpose |
|---|---|
| `01_transcriptomics.R` | RMA normalisation, limma with batch as covariate, QC plots, volcano plot, heatmaps, GO/KEGG enrichment |
| `02_mendelian_randomisation.R` | MR analyses A and B, sensitivity analyses (HLA exclusion, heterogeneity, MR-Egger intercept, Steiger), colocalisation with prior sensitivity |
| `03_forest_plot_MR_A.R` | Forest plot for MR analysis A |
| `04_power_MDE.R` | Post hoc power: minimum detectable effects for analysis A (Table 5) and analysis B |
| `results/DEG/` | Output of script 01: DEG tables, enrichment tables, figures, session information |
| `results/MR/` | Output of scripts 02 to 04: instruments, MR estimates, sensitivity analyses, colocalisation, session information |

## Data

No individual-level data are included. All inputs are public.

- **GSE40611 raw data.** Download `GSE40611_RAW.tar` from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40611), extract it, and place the 49 `*.CEL.gz` files in a folder named `data/` inside the repository. Sample phenotypes (disease status, processing batch) are downloaded automatically from GEO by `GEOquery`.
- **GWAS and eQTL summary statistics.** Retrieved through the [OpenGWAS API](https://api.opengwas.io/) at run time:

| Trait | OpenGWAS ID | Source |
|---|---|---|
| Sjögren's syndrome | `ebi-a-GCST90013879` | UK Biobank (Mbatchou et al., 2021) |
| LDL, HDL, triglycerides | `ieu-a-300`, `ieu-a-299`, `ieu-a-302` | GLGC (Willer et al., 2013) |
| C-reactive protein | `ieu-b-35` | CHARGE (Ligthart et al., 2018) |
| Coronary heart disease | `ieu-a-7` | CARDIoGRAMplusC4D (Nikpay et al., 2015) |
| Systolic, diastolic blood pressure | `ieu-b-38`, `ieu-b-39` | ICBP (Evangelou et al., 2018) |
| Gene expression (blood) | `eqtl-a-<Ensembl ID>` | eQTLGen (Võsa et al., 2021) |

## Requirements

R 4.5.0. Main package versions used: affy 1.86.0, limma 3.64.3, GEOquery 2.76.0, hgu133plus2.db 3.13.0, org.Hs.eg.db 3.21.0, clusterProfiler 4.16.0, TwoSampleMR 0.7.14, ieugwasr 1.2.1, coloc 5.2.3. Full details are in `results/DEG/sessionInfo.txt` and `results/MR/sessionInfo.txt`.

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("affy", "limma", "GEOquery", "hgu133plus2.db", "hgu133plus2cdf",
                       "org.Hs.eg.db", "clusterProfiler", "enrichplot"), update = FALSE)
install.packages(c("pheatmap", "ggrepel", "ggplot2", "dplyr", "tidyr", "patchwork",
                   "ieugwasr", "coloc", "psych"))
remotes::install_github("MRCIEU/TwoSampleMR")
```

### OpenGWAS token

The MR script needs a personal OpenGWAS token. Request one at <https://api.opengwas.io/>, then add it to your `~/.Renviron` file:

```
OPENGWAS_JWT=<your token>
```

Restart R and check the token with `ieugwasr::user()`. Never commit this file or the token to the repository.

## How to run

Open `sjogren-analysis.Rproj` in RStudio (or set the working directory to the repository folder), then run:

```r
source("01_transcriptomics.R")          # -> results/DEG/
source("02_mendelian_randomisation.R")  # -> results/MR/
source("03_forest_plot_MR_A.R")         # -> results/MR/A_forest.png
source("04_power_MDE.R")                # -> results/MR/Power_MDE_A.csv, Power_MDE_B.csv
```

Script 01 needs internet access only to download the GEO phenotype table and for the KEGG step. Script 02 needs internet access and a valid OpenGWAS token. Scripts 03 and 04 read the CSV files written by script 02.

## Reproducibility notes

- The MR results in the paper were generated on 19 November 2026. OpenGWAS data and the LD reference panel may be updated, so instrument lists or estimates can differ slightly in later runs.
- Batch correction (`limma::removeBatchEffect`) is used for PCA and heatmaps only. Differential expression is tested with batch as a covariate in the limma model.

## Licence

Code is released under the MIT Licence (see `LICENSE`).

## Contact

Muhammad Agil Wijaya Faradis, Faculty of Medicine, Universitas Brawijaya, Malang, Indonesia (agilwijaya@student.ub.ac.id)
