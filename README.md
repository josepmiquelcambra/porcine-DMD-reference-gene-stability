# Reference gene stability analysis in a porcine DMD model

This repository contains the R code used to reproduce the reference-gene stability analyses performed in a porcine Duchenne muscular dystrophy (DMD) model.

The workflow evaluates candidate reference genes separately in **skeletal muscle** and **cardiac muscle** and includes quality-control procedures, visualization, individual stability methods, a custom Weighted Stability Index (WSI), a Consensus Stability Score (CSS), and CSS-guided pairwise variation analysis.

> **Reproducibility note:** The analytical implementation in this repository is preserved as used for the reported calculations. The purpose of this repository is to reproduce the analysis, not to retrospectively modify the mathematical procedures.

## Repository contents

- `File_S1_reference_gene_stability_analysis.R` — complete R analysis script.
- `Data_table.xlsx` — raw qPCR Cq measurements required as input.
- `Metadata.xlsx` — sample metadata required as input.
- `Results.xlsx` — generated output containing QC and stability-analysis tables.
- `sessionInfo.txt` — generated record of the R environment and package versions.
- `LICENSE` — repository license.

`Results.xlsx` and `sessionInfo.txt` are created automatically when the script is run.

## Requirements

The analysis requires R and the following packages.

### CRAN packages

- `readxl`
- `dplyr`
- `ggplot2`
- `stringr`
- `tidyr`
- `writexl`

They can be installed with:

```r
install.packages(c(
  "readxl",
  "dplyr",
  "ggplot2",
  "stringr",
  "tidyr",
  "writexl"
))
```

### Bioconductor package

- `NormqPCR`

Installation:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("NormqPCR")
```

The script checks that all required packages are available before starting the analysis.

## Input files

Both input files must be placed in the same working directory as the R script and must retain the filenames shown below.

### `Data_table.xlsx`

This file contains the raw qPCR measurements and requires four columns:

| Column | Description |
|---|---|
| `Sample_ID` | Biological sample identifier |
| `Gene` | Candidate reference gene |
| `Cq` | qPCR Cq value |
| `Replicate` | Technical replicate identifier |

The expected technical-replicate labels are:

- `Rep1`
- `Rep2`

### `Metadata.xlsx`

This file contains sample-level experimental information and requires four columns:

| Column | Description |
|---|---|
| `Sample_ID` | Biological sample identifier |
| `Tissue` | Tissue of sample collection |
| `Genotype` | Experimental genotype |
| `Time` | Sampling/developmental time point |

The analysis expects the following values:

**Tissue**
- `Skeletal Muscle`
- `Cardiac Muscle`

**Genotype**
- `WT`
- `DMD`

**Time**
- `Fetal`
- `1mo`
- `3.5mo`

`Sample_ID` values must be unique in `Metadata.xlsx`, and every sample present in `Data_table.xlsx` must have a corresponding metadata entry.

## Analysis workflow

The script performs the following steps.

### 1. Input validation and preprocessing

The input tables are checked for the required columns, duplicated metadata identifiers, unexpected genotype/time/replicate labels, and unmatched sample identifiers. Text variables are standardized before the qPCR table and metadata are merged.

### 2. Raw Cq visualization

Raw Cq distributions are visualized separately for skeletal and cardiac muscle using genotype-specific boxplots with technical replicate and time-point information.

### 3. Technical replicate variability

For each sample-gene combination, the standard deviation of the technical replicate Cq values is calculated. Mean technical variability is then summarized by gene and tissue.

### 4. Technical-replicate averaging

Technical replicates are averaged for each biological sample and candidate gene:

- `Mean_Cq` = mean Cq of the technical replicates.
- `SD_Cq` = standard deviation of the technical replicates.

### 5. Principal component analysis

PCA is performed separately for skeletal and cardiac muscle using mean Cq values. Variables are centered and scaled before PCA.

### 6. Outlier identification and removal

Two quality-control criteria are applied exactly as used in the study.

**Mean Cq outliers**

A sample-gene observation is flagged when both conditions are met:

- `Mean_Cq` is outside `Q1 - 7 × IQR` or `Q3 + 7 × IQR`, and
- its absolute difference from the tissue/gene median is at least 7 Cq.

**Technical SD outliers**

A sample-gene observation is flagged when:

- `SD_Cq > Q3 + 1.5 × IQR`.

Flagged `Sample_ID + Tissue + Gene` combinations are removed before the stability analyses.

### 7. Relative quantities

For methods requiring relative quantities, the script calculates:

```text
Quantity = 2^-(Mean_Cq - minimum Mean_Cq for that tissue and gene)
```

### 8. Reference-gene stability analyses

The following approaches are applied separately to skeletal and cardiac muscle.

#### BestKeeper-derived analysis

The implementation used in the study reports Cq-based descriptive statistics, identifies genes according to the specified SD criterion, constructs the study-specific BestKeeper-derived index, and calculates Pearson correlations.

The ranking used for downstream consensus analysis is based on `SD_Cq`, with lower values indicating greater stability.

#### geNorm

An iterative geNorm-style procedure calculates gene stability from the mean standard deviation of pairwise log2 expression ratios. Genes are progressively eliminated according to their M values.

#### NormFinder

NormFinder is implemented using:

```r
NormqPCR::stabMeasureRho()
```

Experimental groups are defined by the combination of genotype and time point (`Genotype × Time`).

#### Comparative Delta Ct

For each candidate gene, pairwise Delta Ct values relative to every other candidate reference gene are calculated. The mean of the corresponding standard deviations is used as the stability score.

#### Weighted Stability Index (WSI)

The custom WSI integrates four characteristics of candidate-gene expression:

- standard deviation of mean Cq;
- coefficient of variation of mean Cq;
- Cq range;
- absolute difference in mean Cq between genotypes.

The genotype-difference component contributes 50% of the final WSI score. The SD, CV, and Cq-range components collectively contribute the remaining 50%.

Lower WSI scores indicate greater stability.

#### Consensus Stability Score (CSS)

The quantitative outputs used from the five stability approaches are standardized as z-scores:

- BestKeeper-derived value: `SD_Cq`
- geNorm value: `geNorm_M`
- NormFinder value: `NormFinder_score`
- Comparative Delta Ct value: `deltaCt_score`
- WSI value: `WSI_score`

The CSS is calculated as the arithmetic mean of these standardized values.

Lower CSS values indicate greater overall stability.

#### CSS-guided pairwise variation (Consensus-geNorm V)

Pairwise variation is calculated using candidate genes ordered according to the CSS ranking. Therefore, genes are introduced according to the consensus ranking rather than the native geNorm ranking.

The resulting `Vn/n+1` values evaluate the effect of adding the next CSS-ranked candidate reference gene to the normalization factor.

## Output

Running the script creates `Results.xlsx`, containing the following worksheets:

### Quality-control tables

- `Outlier_summary_by_sample`
- `Mean_Cq_outliers`
- `Technical_SD_outliers`

### Stability analyses

- `BestKeeper_SkeletalMuscle`
- `BestKeeper_CardiacMuscle`
- `geNorm_SkeletalMuscle`
- `geNorm_CardiacMuscle`
- `NormFinder_SkeletalMuscle`
- `NormFinder_CardiacMuscle`
- `DeltaCt_SkeletalMuscle`
- `DeltaCt_CardiacMuscle`
- `WSI_SkeletalMuscle`
- `WSI_CardiacMuscle`
- `CSS_SkeletalMuscle`
- `CSS_CardiacMuscle`

### Pairwise variation

- `Consensus_geNorm_V_SkeletalMuscle`
- `Consensus_geNorm_V_CardiacMuscle`

The script also creates `sessionInfo.txt`, which records the R version, operating system, and loaded package versions used for the analysis.

The plots generated during the workflow are displayed in the active R graphics device but are not automatically exported as separate files.

## Running the analysis

1. Download or clone this repository.
2. Place `Data_table.xlsx` and `Metadata.xlsx` in the same working directory as `File_S1_reference_gene_stability_analysis.R`.
3. Install the required R packages.
4. Set the repository folder as the R working directory.
5. Run the complete script.

For example:

```r
setwd("path/to/porcine-DMD-reference-gene-stability")
source("File_S1_reference_gene_stability_analysis.R")
```

After successful completion, the working directory will contain:

```text
Results.xlsx
sessionInfo.txt
```

## Reproducibility

This repository is intended to provide the exact analysis workflow associated with the study. Methodological choices and calculations in the analysis script are intentionally retained as used for the reported results.

For long-term reproducibility, the version of the repository associated with the publication should be archived as a fixed release and assigned a DOI through Zenodo.

## Citation

If you use this code, please cite the associated publication and the archived Zenodo software release.

**Zenodo DOI:** https://doi.org/10.5281/zenodo.22706354

## License

This repository is distributed under the terms of the `LICENSE` file included in the repository.
