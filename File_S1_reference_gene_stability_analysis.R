# ============================================================================
# File S1. Reference gene stability analysis in a porcine DMD model
# ============================================================================
#
# Purpose
# -------
# This script reproduces the reference-gene stability analyses used in the
# associated study.
#
# Required input files (placed in the working directory)
# ------------------------------------------------------
#   1. Data_table.xlsx
#      Required columns:
#        - Sample_ID : biological sample identifier
#        - Gene      : candidate reference gene
#        - Cq        : qPCR Cq value
#        - Replicate : technical replicate identifier (Rep1 or Rep2)
#
#   2. Metadata.xlsx
#      Required columns:
#        - Sample_ID : biological sample identifier
#        - Tissue    : tissue of sample collection
#        - Genotype  : WT or DMD
#        - Time      : Fetal, 1mo, or 3.5mo
#
# Main output
# -----------
#   Results.xlsx      : QC tables and all algorithm-specific results
#   sessionInfo.txt   : R and package versions used to run the script
#
# Analyses included
# -----------------
#   - Raw Cq visualization
#   - Technical replicate variability
#   - PCA
#   - Outlier identification and removal
#   - BestKeeper-derived analysis
#   - geNorm
#   - NormFinder
#   - Comparative Delta Ct
#   - Weighted Stability Index (WSI)
#   - Consensus Stability Score (CSS)
#   - CSS-guided pairwise variation (Consensus-geNorm V)

# ============================================================
# Required packages
# ============================================================

required_cran_packages <- c(
  "readxl",
  "dplyr",
  "ggplot2",
  "stringr",
  "tidyr",
  "writexl"
)

missing_cran_packages <- required_cran_packages[
  !vapply(required_cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran_packages) > 0) {
  stop(
    "Missing required CRAN package(s): ",
    paste(missing_cran_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

if (!requireNamespace("NormqPCR", quietly = TRUE)) {
  stop(
    "Missing required Bioconductor package: NormqPCR. ",
    "Install it before running this script, for example with ",
    "BiocManager::install('NormqPCR')."
  )
}

library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(writexl)
library(NormqPCR)

# ============================================================
# Read and validate input files
# ============================================================

input_data_file <- "Data_table.xlsx"
input_metadata_file <- "Metadata.xlsx"

if (!file.exists(input_data_file)) {
  stop("Input file not found: ", input_data_file)
}

if (!file.exists(input_metadata_file)) {
  stop("Input file not found: ", input_metadata_file)
}

# Data_table.xlsx contains the raw qPCR measurements.
data_raw <- read_excel(input_data_file)

# Metadata.xlsx contains the sample-level experimental information.
metadata <- read_excel(input_metadata_file)

required_data_columns <- c("Sample_ID", "Gene", "Cq", "Replicate")
required_metadata_columns <- c("Sample_ID", "Tissue", "Genotype", "Time")

missing_data_columns <- setdiff(required_data_columns, colnames(data_raw))
missing_metadata_columns <- setdiff(required_metadata_columns, colnames(metadata))

if (length(missing_data_columns) > 0) {
  stop(
    "Data_table.xlsx is missing required column(s): ",
    paste(missing_data_columns, collapse = ", ")
  )
}

if (length(missing_metadata_columns) > 0) {
  stop(
    "Metadata.xlsx is missing required column(s): ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

# Clean text variables before merging
data_raw <- data_raw %>%
  mutate(
    Sample_ID = str_squish(as.character(Sample_ID)),
    Gene = str_squish(as.character(Gene)),
    Replicate = str_squish(as.character(Replicate)),
    Cq = as.numeric(Cq)
  )

metadata <- metadata %>%
  mutate(
    Sample_ID = str_squish(as.character(Sample_ID)),
    Tissue = str_squish(as.character(Tissue)),
    Genotype = str_squish(as.character(Genotype)),
    Time = str_squish(as.character(Time))
  )

# Validate identifiers and categorical values before merging/factor conversion.
# These checks do not modify the data; they prevent silent execution with an
# input structure that differs from that used in the study.
if (anyDuplicated(metadata$Sample_ID)) {
  stop("Duplicated Sample_ID values detected in Metadata.xlsx.")
}

expected_genotypes <- c("WT", "DMD")
expected_times <- c("Fetal", "1mo", "3.5mo")
expected_replicates <- c("Rep1", "Rep2")

unexpected_genotypes <- setdiff(unique(metadata$Genotype), expected_genotypes)
unexpected_times <- setdiff(unique(metadata$Time), expected_times)
unexpected_replicates <- setdiff(unique(data_raw$Replicate), expected_replicates)

if (length(unexpected_genotypes) > 0) {
  stop(
    "Unexpected Genotype value(s): ",
    paste(unexpected_genotypes, collapse = ", ")
  )
}

if (length(unexpected_times) > 0) {
  stop(
    "Unexpected Time value(s): ",
    paste(unexpected_times, collapse = ", ")
  )
}

if (length(unexpected_replicates) > 0) {
  stop(
    "Unexpected Replicate value(s): ",
    paste(unexpected_replicates, collapse = ", ")
  )
}

# Merge raw qPCR data with sample metadata.
data_merged <- data_raw %>%
  left_join(metadata, by = "Sample_ID")

# Check whether every qPCR observation was matched to metadata.
unmatched_samples <- data_merged %>%
  filter(is.na(Tissue) | is.na(Genotype) | is.na(Time))

if (nrow(unmatched_samples) > 0) {
  stop(
    "One or more Sample_ID values in Data_table.xlsx could not be matched ",
    "to Metadata.xlsx."
  )
}

# Define factor order exactly as used in the study.
data_merged <- data_merged %>%
  mutate(
    Genotype = factor(Genotype, levels = c("WT", "DMD")),
    Time = factor(Time, levels = c("Fetal", "1mo", "3.5mo")),
    Replicate = factor(Replicate, levels = c("Rep1", "Rep2"))
  )

# ============================================================
# BOXPLOT Cq values Skeletal Muscle
# ============================================================

boxplotCq_skeletalmuscle <- data_merged %>%
  filter(Tissue == "Skeletal Muscle") %>%
  ggplot(aes(x = Gene, y = Cq)) +
  geom_boxplot(
    aes(
      color = Genotype,
      group = interaction(Gene, Genotype)
    ),
    position = position_dodge(width = 0.50),
    width = 0.30,
    fill = NA,
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  geom_point(
    aes(
      color = Genotype,
      fill = Replicate,
      shape = Time,
      group = Genotype
    ),
    position = position_jitterdodge(
      jitter.width = 0.08,
      jitter.height = 0,
      dodge.width = 0.50
    ),
    size = 1.6,
    alpha = 0.85,
    stroke = 0.7
  ) +
  scale_color_manual(values = c("WT" = "blue", "DMD" = "red")) +
  scale_fill_manual(values = c("Rep1" = "black", "Rep2" = "white")) +
  scale_shape_manual(values = c("Fetal" = 21, "1mo" = 24, "3.5mo" = 22)) +
  labs(
    title = "Figure 1A. Raw Cq values in skeletal muscle",
    x = "Candidate housekeeping gene",
    y = "Cq value",
    color = "Group",
    fill = "Technical replicate",
    shape = "Time point"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

boxplotCq_skeletalmuscle

# ============================================================
# BOXPLOT Cq values Cardiac Muscle
# ============================================================

boxplotCq_cardiacmuscle <- data_merged %>%
  filter(Tissue == "Cardiac Muscle") %>%
  ggplot(aes(x = Gene, y = Cq)) +
  geom_boxplot(
    aes(
      color = Genotype,
      group = interaction(Gene, Genotype)
    ),
    position = position_dodge(width = 0.50),
    width = 0.30,
    fill = NA,
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  geom_point(
    aes(
      color = Genotype,
      fill = Replicate,
      shape = Time,
      group = Genotype
    ),
    position = position_jitterdodge(
      jitter.width = 0.08,
      jitter.height = 0,
      dodge.width = 0.50
    ),
    size = 1.6,
    alpha = 0.85,
    stroke = 0.7
  ) +
  scale_color_manual(values = c("WT" = "blue", "DMD" = "red")) +
  scale_fill_manual(values = c("Rep1" = "black", "Rep2" = "white")) +
  scale_shape_manual(values = c("Fetal" = 21, "1mo" = 24, "3.5mo" = 22)) +
  labs(
    title = "Figure 1B. Raw Cq values in cardiac muscle",
    x = "Candidate housekeeping gene",
    y = "Cq value",
    color = "Group",
    fill = "Technical replicate",
    shape = "Time point"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

boxplotCq_cardiacmuscle

# ============================================================
# Technical deviation analyses
# ============================================================

tech_sd <- data_merged %>%
  group_by(Sample_ID, Gene, Tissue) %>%
  summarise(
    SD = sd(Cq, na.rm = TRUE),
    .groups = "drop"
  )

tech_sd_summary <- tech_sd %>%
  group_by(Gene, Tissue) %>%
  summarise(
    mean_SD = mean(SD, na.rm = TRUE),
    .groups = "drop"
  )

SD_skeletalmuscle <- tech_sd_summary %>%
  filter(Tissue == "Skeletal Muscle") %>%
  ggplot(aes(x = Gene, y = mean_SD)) +
  geom_bar(stat = "identity", fill = "grey60") +
  labs(
    title = "SD_skeletalmuscle. Technical variability (SD) in skeletal muscle",
    x = "Candidate housekeeping gene",
    y = "Mean SD (Cq)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

SD_cardiacmuscle <- tech_sd_summary %>%
  filter(Tissue == "Cardiac Muscle") %>%
  ggplot(aes(x = Gene, y = mean_SD)) +
  geom_bar(stat = "identity", fill = "grey60") +
  labs(
    title = "SD_cardiacmuscle. Technical variability (SD) in cardiac muscle",
    x = "Candidate housekeeping gene",
    y = "Mean SD (Cq)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

SD_skeletalmuscle
SD_cardiacmuscle

# ============================================================
# Average technical replicates
# ============================================================
# Mean Cq and technical-replicate SD are calculated for each biological
# sample and candidate gene exactly as used in the original analysis.

data_mean <- data_merged %>%
  dplyr::group_by(Sample_ID, Tissue, Genotype, Time, Gene) %>%
  dplyr::summarise(
    Mean_Cq = mean(Cq, na.rm = TRUE),
    SD_Cq = sd(Cq, na.rm = TRUE),
    .groups = "drop"
  )
# ============================================================
# PCA analyses
# ============================================================

calculate_pca_plot <- function(data_tissue, tissue_name) {
  
  pca_data <- data_tissue %>%
    dplyr::select(Sample_ID, Genotype, Time, Gene, Mean_Cq) %>%
    tidyr::pivot_wider(
      names_from = Gene,
      values_from = Mean_Cq
    )
  
  metadata <- pca_data %>%
    dplyr::select(Sample_ID, Genotype, Time)
  
  expression_matrix <- pca_data %>%
    dplyr::select(-Sample_ID, -Genotype, -Time)
  
  pca <- prcomp(
    expression_matrix,
    center = TRUE,
    scale. = TRUE
  )
  
  variance <- summary(pca)$importance[2,] * 100
  
  scores <- as.data.frame(pca$x)
  
  scores <- cbind(metadata, scores)
  
  p <- ggplot(
    scores,
    aes(
      PC1,
      PC2,
      colour = Genotype,
      shape = Time
    )
  ) +
    
    geom_point(size = 3.5) +
    
    stat_ellipse(
      aes(group = Genotype),
      linewidth = 0.8,
      alpha = 0.3
    ) +
    
    scale_colour_manual(
      values = c(
        WT = "blue",
        DMD = "red"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        Fetal = 16,
        `1mo` = 17,
        `3.5mo` = 15
      )
    ) +
    
    labs(
      title = paste("PCA -", tissue_name),
      x = paste0("PC1 (", round(variance[1],1), "%)"),
      y = paste0("PC2 (", round(variance[2],1), "%)")
    ) +
    
    theme_classic() +
    
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
  
  return(p)
}

pca_skeletalmuscle <- calculate_pca_plot(
  data_mean %>% filter(Tissue == "Skeletal Muscle"),
  "Skeletal Muscle"
)

pca_cardiacmuscle <- calculate_pca_plot(
  data_mean %>% filter(Tissue == "Cardiac Muscle"),
  "Cardiac Muscle"
)

pca_skeletalmuscle
pca_cardiacmuscle

# ============================================================
# OUTLIER IDENTIFICATION
# ============================================================
# Mean Cq outlier criterion used in the study:
#   - Mean_Cq outside Q1 - 7*IQR or Q3 + 7*IQR, AND
#   - absolute difference from the tissue/gene median >= 7 Cq.
#
# Technical SD outlier criterion used in the study:
#   - SD_Cq > Q3 + 1.5*IQR.
#
# These thresholds are intentionally preserved exactly as used for the
# reported calculations.

outliers_cq <- data_mean %>%
  group_by(Tissue, Gene) %>%
  mutate(
    median_Cq = median(Mean_Cq, na.rm = TRUE),
    Q1 = quantile(Mean_Cq, 0.25, na.rm = TRUE),
    Q3 = quantile(Mean_Cq, 0.75, na.rm = TRUE),
    IQR_value = IQR(Mean_Cq, na.rm = TRUE),
    lower_limit = Q1 - 7 * IQR_value,
    upper_limit = Q3 + 7 * IQR_value,
    abs_diff_median = abs(Mean_Cq - median_Cq),
    Cq_outlier =
      (Mean_Cq < lower_limit | Mean_Cq > upper_limit) &
      abs_diff_median >= 7
  ) %>%
  ungroup() %>%
  filter(Cq_outlier == TRUE) %>%
  dplyr::select(
    Sample_ID, Tissue, Genotype, Time, Gene,
    Mean_Cq, SD_Cq,
    median_Cq, abs_diff_median,
    lower_limit, upper_limit
  )

outliers_technical_sd <- data_mean %>%
  filter(!is.na(SD_Cq)) %>%
  group_by(Tissue, Gene) %>%
  mutate(
    Q1_SD = quantile(SD_Cq, 0.25, na.rm = TRUE),
    Q3_SD = quantile(SD_Cq, 0.75, na.rm = TRUE),
    IQR_SD = IQR(SD_Cq, na.rm = TRUE),
    upper_limit_SD = Q3_SD + 1.5 * IQR_SD,
    SD_outlier = SD_Cq > upper_limit_SD
  ) %>%
  ungroup() %>%
  filter(SD_outlier == TRUE) %>%
  dplyr::select(
    Sample_ID, Tissue, Genotype, Time, Gene,
    Mean_Cq, SD_Cq,
    upper_limit_SD
  )

outlier_summary_by_sample <- bind_rows(
  outliers_cq %>%
    mutate(Outlier_type = "Mean_Cq_outlier") %>%
    dplyr::select(Sample_ID, Tissue, Genotype, Time, Gene, Outlier_type),
  
  outliers_technical_sd %>%
    mutate(Outlier_type = "Technical_SD_outlier") %>%
    dplyr::select(Sample_ID, Tissue, Genotype, Time, Gene, Outlier_type)
) %>%
  group_by(Sample_ID, Tissue, Genotype, Time, Outlier_type) %>%
  summarise(
    n_flagged_genes = n_distinct(Gene),
    flagged_genes = paste(unique(Gene), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(Tissue, Sample_ID, Outlier_type)



# ============================================================
# REMOVE OUTLIERS BEFORE ALGORITHM ANALYSES
# Removes flagged Sample_ID + Tissue + Gene combinations
# from both Mean_Cq outliers and Technical_SD outliers
# ============================================================

outliers_to_remove <- dplyr::bind_rows(
  outliers_cq %>%
    dplyr::select(Sample_ID, Tissue, Gene),
  
  outliers_technical_sd %>%
    dplyr::select(Sample_ID, Tissue, Gene)
) %>%
  dplyr::distinct()

data_mean_clean <- data_mean %>%
  dplyr::anti_join(
    outliers_to_remove,
    by = c("Sample_ID", "Tissue", "Gene")
  )

data_mean_q_clean <- data_mean_clean %>%
  dplyr::group_by(Tissue, Gene) %>%
  dplyr::mutate(
    min_Cq = min(Mean_Cq, na.rm = TRUE),
    Quantity = 2^-(Mean_Cq - min_Cq)
  ) %>%
  dplyr::ungroup()

# ============================================================
# ALGORITHM ANALYSES
# ============================================================
# 1. BestKeeper-derived analysis = Cq-based descriptive statistics plus the
#    study-specific BestKeeper index implementation documented below.
# 2. geNorm = relative quantities (RQmin).
# 3. NormFinder = relative quantities (RQmin), grouped by Genotype x Time.
# 4. Comparative Delta Ct = Cq.
# 5. Weighted Stability Index (WSI).
# 6. Consensus Stability Score (CSS).
# ============================================================

# ============================================================
# 1. BESTKEEPER
# ============================================================
# Implementation preserved exactly as used in the study:
#   - Descriptive statistics are calculated from Mean_Cq values.
#   - Genes with SD_Cq <= 1 are included in the index.
#   - If fewer than two genes meet this threshold, the three genes with the
#     lowest SD_Cq (or all genes if fewer than three are available) are used.
#   - For index construction, Cq values are transformed as 2^-Cq and the
#     geometric mean is calculated across included genes.
#   - Pearson correlations are calculated between each gene's 2^-Cq values
#     and this index.
#   - BestKeeper_rank is ordered by SD_Cq.

calculate_bestkeeper <- function(data_tissue, sd_cutoff = 1) {
  
  cq_wide <- data_tissue %>%
    dplyr::select(Sample_ID, Gene, Mean_Cq) %>%
    tidyr::pivot_wider(names_from = Gene, values_from = Mean_Cq)
  
  cq_matrix <- cq_wide %>%
    dplyr::select(-Sample_ID) %>%
    as.data.frame()
  
  genes <- colnames(cq_matrix)
  
  initial_stats <- data.frame(
    Gene = genes,
    mean_Cq = sapply(cq_matrix, mean, na.rm = TRUE),
    SD_Cq = sapply(cq_matrix, sd, na.rm = TRUE),
    CV_percent = sapply(cq_matrix, function(x) {
      sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
    })
  )
  
  stable_genes <- initial_stats %>%
    dplyr::filter(SD_Cq <= sd_cutoff) %>%
    dplyr::pull(Gene)
  
  if (length(stable_genes) < 2) {
    stable_genes <- initial_stats %>%
      dplyr::arrange(SD_Cq) %>%
      dplyr::slice_head(n = min(3, n())) %>%
      dplyr::pull(Gene)
  }
  
  expr_matrix <- 2^(-cq_matrix[, stable_genes, drop = FALSE])
  
  bestkeeper_index <- apply(expr_matrix, 1, function(x) {
    exp(mean(log(x), na.rm = TRUE))
  })
  
  results <- initial_stats %>%
    dplyr::mutate(
      Included_in_index = Gene %in% stable_genes,
      Pearson_r = sapply(Gene, function(gene) {
        gene_expr <- 2^(-cq_matrix[[gene]])
        cor(
          gene_expr,
          bestkeeper_index,
          use = "pairwise.complete.obs",
          method = "pearson"
        )
      }),
      Pearson_p = sapply(Gene, function(gene) {
        gene_expr <- 2^(-cq_matrix[[gene]])
        complete_cases <- complete.cases(gene_expr, bestkeeper_index)
        cor.test(
          gene_expr[complete_cases],
          bestkeeper_index[complete_cases],
          method = "pearson"
        )$p.value
      })
    ) %>%
    dplyr::arrange(SD_Cq) %>%
    dplyr::mutate(BestKeeper_rank = dplyr::row_number())
  
  return(results)
}

bestkeeper_skeletalmuscle <- data_mean_clean %>%
  dplyr::filter(Tissue == "Skeletal Muscle") %>%
  calculate_bestkeeper(sd_cutoff = 1)

bestkeeper_cardiacmuscle <- data_mean_clean %>%
  dplyr::filter(Tissue == "Cardiac Muscle") %>%
  calculate_bestkeeper(sd_cutoff = 1)

# ============================================================
# 2. geNorm
# ============================================================
# Iterative elimination is based on the mean SD of pairwise log2 expression
# ratios. The implementation below is preserved exactly as used in the study.

calculate_genorm_iterative <- function(data_tissue) {
  
  q_wide <- data_tissue %>%
    dplyr::select(Sample_ID, Gene, Quantity) %>%
    tidyr::pivot_wider(names_from = Gene, values_from = Quantity)
  
  q_matrix <- q_wide %>%
    dplyr::select(-Sample_ID) %>%
    as.data.frame()
  
  remaining_genes <- colnames(q_matrix)
  ranking <- data.frame()
  
  while (length(remaining_genes) > 2) {
    
    M_values <- sapply(remaining_genes, function(gene) {
      other_genes <- remaining_genes[remaining_genes != gene]
      
      pairwise_sd <- sapply(other_genes, function(other_gene) {
        ratio <- log2(q_matrix[[gene]] / q_matrix[[other_gene]])
        sd(ratio, na.rm = TRUE)
      })
      
      mean(pairwise_sd, na.rm = TRUE)
    })
    
    worst_gene <- names(which.max(M_values))
    
    ranking <- rbind(
      ranking,
      data.frame(
        Gene = worst_gene,
        geNorm_M = M_values[worst_gene],
        elimination_order = length(remaining_genes)
      )
    )
    
    remaining_genes <- remaining_genes[remaining_genes != worst_gene]
  }
  
  final_M <- sapply(remaining_genes, function(gene) {
    other_gene <- remaining_genes[remaining_genes != gene]
    ratio <- log2(q_matrix[[gene]] / q_matrix[[other_gene]])
    sd(ratio, na.rm = TRUE)
  })
  
  ranking <- rbind(
    ranking,
    data.frame(
      Gene = remaining_genes,
      geNorm_M = final_M,
      elimination_order = 2
    )
  )
  
  ranking <- ranking %>%
    dplyr::arrange(elimination_order) %>%
    dplyr::mutate(geNorm_rank = dplyr::row_number())
  
  return(ranking)
}

genorm_skeletalmuscle <- data_mean_q_clean %>%
  dplyr::filter(Tissue == "Skeletal Muscle") %>%
  calculate_genorm_iterative()

genorm_cardiacmuscle <- data_mean_q_clean %>%
  dplyr::filter(Tissue == "Cardiac Muscle") %>%
  calculate_genorm_iterative()


# ============================================================
# 3. NormFinder
# ============================================================
# NormFinder is implemented with NormqPCR::stabMeasureRho().
# Rows of rq_matrix correspond to biological samples and columns to candidate
# reference genes. Experimental groups are defined as Genotype x Time.

calculate_normfinder <- function(data_tissue) {
  
  rq_wide <- data_tissue %>%
    dplyr::select(Sample_ID, Genotype, Time, Gene, Quantity) %>%
    tidyr::pivot_wider(names_from = Gene, values_from = Quantity)
  
  sample_info <- rq_wide %>%
    dplyr::select(Sample_ID, Genotype, Time) %>%
    dplyr::mutate(
      NormFinder_group = paste(Genotype, Time, sep = "_")
    )
  
  rq_matrix <- rq_wide %>%
    dplyr::select(-Sample_ID, -Genotype, -Time) %>%
    as.matrix()
  
  rownames(rq_matrix) <- as.character(rq_wide$Sample_ID)
  
  group_factor <- as.factor(sample_info$NormFinder_group)
  
  nf <- NormqPCR::stabMeasureRho(
    rq_matrix,
    group = group_factor,
    log = FALSE,
    na.rm = TRUE,
    returnAll = TRUE
  )
  
  normfinder_results <- data.frame(
    Gene = colnames(rq_matrix),
    NormFinder_score = as.numeric(nf$rho)
  ) %>%
    dplyr::arrange(NormFinder_score) %>%
    dplyr::mutate(NormFinder_rank = dplyr::row_number())
  
  return(normfinder_results)
}

normfinder_skeletalmuscle <- data_mean_q_clean %>%
  dplyr::filter(Tissue == "Skeletal Muscle") %>%
  calculate_normfinder()

normfinder_cardiacmuscle <- data_mean_q_clean %>%
  dplyr::filter(Tissue == "Cardiac Muscle") %>%
  calculate_normfinder()

# ============================================================
# 4. Comparative Delta Ct
# ============================================================
# For each candidate gene, the SD of its Delta Ct relative to every other
# candidate gene is calculated and then averaged, exactly as used in the study.

calculate_deltaCt <- function(data_tissue) {
  
  cq_wide <- data_tissue %>%
    dplyr::select(Sample_ID, Gene, Mean_Cq) %>%
    tidyr::pivot_wider(names_from = Gene, values_from = Mean_Cq)
  
  cq_matrix <- cq_wide %>%
    dplyr::select(-Sample_ID) %>%
    as.data.frame()
  
  genes <- colnames(cq_matrix)
  
  results <- data.frame(
    Gene = genes,
    deltaCt_score = NA_real_
  )
  
  for (gene in genes) {
    
    other_genes <- genes[genes != gene]
    
    pairwise_sd <- sapply(other_genes, function(other_gene) {
      delta_ct <- cq_matrix[[gene]] - cq_matrix[[other_gene]]
      sd(delta_ct, na.rm = TRUE)
    })
    
    results$deltaCt_score[results$Gene == gene] <- mean(pairwise_sd, na.rm = TRUE)
  }
  
  results <- results %>%
    dplyr::arrange(deltaCt_score) %>%
    dplyr::mutate(deltaCt_rank = dplyr::row_number())
  
  return(results)
}

deltaCt_skeletalmuscle <- data_mean_clean %>%
  dplyr::filter(Tissue == "Skeletal Muscle") %>%
  calculate_deltaCt()

deltaCt_cardiacmuscle <- data_mean_clean %>%
  dplyr::filter(Tissue == "Cardiac Muscle") %>%
  calculate_deltaCt()


# ============================================================
# 5. Weighted Stability Index (WSI)
# ============================================================
# Custom stability metric used in this study.
# Components:
#   - SD of Mean_Cq
#   - CV of Mean_Cq
#   - Cq range
#   - absolute difference in mean Cq between genotypes
#
# The genotype component contributes 50% of the final WSI score. SD, CV and
# Cq range collectively contribute the remaining 50%. Lower WSI scores indicate
# greater stability. The weighting and rank-based calculation are preserved
# exactly as used in the study.

calculate_wsi <- function(data_tissue, genotype_col = "Genotype") {
  
  global_stats <- data_tissue %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      WSI_SD_Cq = sd(Mean_Cq, na.rm = TRUE),
      WSI_mean_Cq = mean(Mean_Cq, na.rm = TRUE),
      WSI_CV_percent = (WSI_SD_Cq / WSI_mean_Cq) * 100,
      WSI_Cq_range = max(Mean_Cq, na.rm = TRUE) - min(Mean_Cq, na.rm = TRUE),
      .groups = "drop"
    )
  
  genotype_stats <- data_tissue %>%
    dplyr::group_by(Gene, .data[[genotype_col]]) %>%
    dplyr::summarise(
      mean_Cq_genotype = mean(Mean_Cq, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = .data[[genotype_col]],
      values_from = mean_Cq_genotype
    )
  
  genotype_levels <- setdiff(colnames(genotype_stats), "Gene")
  
  if (length(genotype_levels) != 2) {
    stop("WSI method requires exactly two genotype groups.")
  }
  
  genotype_stats <- genotype_stats %>%
    dplyr::mutate(
      WSI_abs_delta_genotype =
        abs(.data[[genotype_levels[1]]] - .data[[genotype_levels[2]]])
    ) %>%
    dplyr::select(Gene, WSI_abs_delta_genotype)
  
  wsi_results <- global_stats %>%
    dplyr::left_join(genotype_stats, by = "Gene") %>%
    dplyr::mutate(
      WSI_SD_rank = rank(WSI_SD_Cq, ties.method = "average"),
      WSI_CV_rank = rank(WSI_CV_percent, ties.method = "average"),
      WSI_range_rank = rank(WSI_Cq_range, ties.method = "average"),
      WSI_genotype_rank = rank(WSI_abs_delta_genotype, ties.method = "average"),
      
      WSI_score =
        ((WSI_SD_rank + WSI_CV_rank + WSI_range_rank) / 3) * 0.5 +
        WSI_genotype_rank * 0.5,
      
      WSI_rank = rank(WSI_score, ties.method = "average")
    ) %>%
    dplyr::arrange(WSI_rank)
  
  return(wsi_results)
}

wsi_skeletalmuscle <- data_mean_clean %>%
  dplyr::filter(Tissue == "Skeletal Muscle") %>%
  calculate_wsi(genotype_col = "Genotype")

wsi_cardiacmuscle <- data_mean_clean %>%
  dplyr::filter(Tissue == "Cardiac Muscle") %>%
  calculate_wsi(genotype_col = "Genotype")


# ============================================================
# 6. Consensus Stability Score (CSS)
# ============================================================
# The quantitative output used from each method is standardized using z-scores:
#   BestKeeper value = SD_Cq
#   geNorm value             = geNorm_M
#   NormFinder value         = NormFinder_score
#   Delta Ct value           = deltaCt_score
#   WSI value                = WSI_score
#
# CSS_score is the arithmetic mean of these standardized values (na.rm = TRUE).
# Lower CSS scores indicate greater overall stability. This calculation is
# preserved exactly as used in the study.

calculate_css_zscore_ranking <- function(
    bestkeeper,
    genorm,
    normfinder,
    deltaCt,
    wsi
) {
  
  ranking_z <- bestkeeper %>%
    dplyr::select(Gene, BestKeeper_value = SD_Cq) %>%
    dplyr::left_join(
      genorm %>% dplyr::select(Gene, geNorm_value = geNorm_M),
      by = "Gene"
    ) %>%
    dplyr::left_join(
      normfinder %>% dplyr::select(Gene, NormFinder_value = NormFinder_score),
      by = "Gene"
    ) %>%
    dplyr::left_join(
      deltaCt %>% dplyr::select(Gene, deltaCt_value = deltaCt_score),
      by = "Gene"
    ) %>%
    dplyr::left_join(
      wsi %>% dplyr::select(Gene, WSI_value = WSI_score),
      by = "Gene"
    ) %>%
    dplyr::mutate(
      BestKeeper_z = as.numeric(scale(BestKeeper_value)),
      geNorm_z = as.numeric(scale(geNorm_value)),
      NormFinder_z = as.numeric(scale(NormFinder_value)),
      deltaCt_z = as.numeric(scale(deltaCt_value)),
      WSI_z = as.numeric(scale(WSI_value))
    ) %>%
    dplyr::mutate(
      CSS_score = rowMeans(
        dplyr::pick(
          BestKeeper_z,
          geNorm_z,
          NormFinder_z,
          deltaCt_z,
          WSI_z
        ),
        na.rm = TRUE
      )
    ) %>%
    dplyr::arrange(CSS_score) %>%
    dplyr::mutate(CSS_rank = dplyr::row_number())
  
  return(ranking_z)
}

css_skeletalmuscle <- calculate_css_zscore_ranking(
  bestkeeper = bestkeeper_skeletalmuscle,
  genorm = genorm_skeletalmuscle,
  normfinder = normfinder_skeletalmuscle,
  deltaCt = deltaCt_skeletalmuscle,
  wsi = wsi_skeletalmuscle
)

css_cardiacmuscle <- calculate_css_zscore_ranking(
  bestkeeper = bestkeeper_cardiacmuscle,
  genorm = genorm_cardiacmuscle,
  normfinder = normfinder_cardiacmuscle,
  deltaCt = deltaCt_cardiacmuscle,
  wsi = wsi_cardiacmuscle
)

# ============================================================
# 7. CSS-guided pairwise variation ("Consensus-geNorm V")
# ============================================================
# Pairwise variation is calculated using genes ordered according to CSS_rank.
# Thus, genes are introduced according to the consensus ranking rather than
# the native geNorm ranking. The implementation below is preserved exactly as
# used in the study.
#
# The general function accepts any ranking table, but in this script it is
# applied to the CSS ranking for both tissues.
# ============================================================

calculate_pairwise_V <- function(data_tissue, ranking_table, rank_column, method_name) {
  
  q_wide <- data_tissue %>%
    dplyr::select(Sample_ID, Gene, Quantity) %>%
    tidyr::pivot_wider(names_from = Gene, values_from = Quantity)
  
  q_matrix <- q_wide %>%
    dplyr::select(-Sample_ID) %>%
    as.data.frame()
  
  ordered_genes <- ranking_table %>%
    dplyr::arrange(.data[[rank_column]]) %>%
    dplyr::pull(Gene)
  
  V_values <- data.frame()
  
  for (i in 2:(length(ordered_genes) - 1)) {
    
    genes_n <- ordered_genes[1:i]
    genes_n1 <- ordered_genes[1:(i + 1)]
    
    NF_n <- apply(q_matrix[, genes_n, drop = FALSE], 1, function(x) {
      exp(mean(log(x), na.rm = TRUE))
    })
    
    NF_n1 <- apply(q_matrix[, genes_n1, drop = FALSE], 1, function(x) {
      exp(mean(log(x), na.rm = TRUE))
    })
    
    V <- sd(log2(NF_n / NF_n1), na.rm = TRUE)
    
    V_values <- rbind(
      V_values,
      data.frame(
        Method = method_name,
        Comparison = paste0("V", i, "/", i + 1),
        n_genes = i,
        Genes_NF_n = paste(genes_n, collapse = " + "),
        Genes_NF_n1 = paste(genes_n1, collapse = " + "),
        V = V
      )
    )
  }
  
  return(V_values)
}


# Consensus-geNorm V: pairwise variation calculated using the CSS ranking.

V_consensus_skeletalmuscle <- data_mean_q_clean %>%
  dplyr::filter(Tissue == "Skeletal Muscle") %>%
  calculate_pairwise_V(
    ranking_table = css_skeletalmuscle,
    rank_column = "CSS_rank",
    method_name = "Consensus_geNorm"
  )

V_consensus_cardiacmuscle <- data_mean_q_clean %>%
  dplyr::filter(Tissue == "Cardiac Muscle") %>%
  calculate_pairwise_V(
    ranking_table = css_cardiacmuscle,
    rank_column = "CSS_rank",
    method_name = "Consensus_geNorm"
  )



# ============================================================
# 8. Export results
# ============================================================

algorithm_tables <- list(
  # Outlier information first in the final Excel file
  "Outlier_summary_by_sample" = outlier_summary_by_sample,
  "Mean_Cq_outliers" = outliers_cq,
  "Technical_SD_outliers" = outliers_technical_sd,
  
  # Algorithm-specific results
  "BestKeeper_SkeletalMuscle" = bestkeeper_skeletalmuscle,
  "BestKeeper_CardiacMuscle" = bestkeeper_cardiacmuscle,
  "geNorm_SkeletalMuscle" = genorm_skeletalmuscle,
  "geNorm_CardiacMuscle" = genorm_cardiacmuscle,
  "NormFinder_SkeletalMuscle" = normfinder_skeletalmuscle,
  "NormFinder_CardiacMuscle" = normfinder_cardiacmuscle,
  "DeltaCt_SkeletalMuscle" = deltaCt_skeletalmuscle,
  "DeltaCt_CardiacMuscle" = deltaCt_cardiacmuscle,
  "WSI_SkeletalMuscle" = wsi_skeletalmuscle,
  "WSI_CardiacMuscle" = wsi_cardiacmuscle,
  "CSS_SkeletalMuscle" = css_skeletalmuscle,
  "CSS_CardiacMuscle" = css_cardiacmuscle,
  
  # Consensus-geNorm V
  "Consensus_geNorm_V_SkeletalMuscle" = V_consensus_skeletalmuscle,
  "Consensus_geNorm_V_CardiacMuscle" = V_consensus_cardiacmuscle
)

writexl::write_xlsx(
  algorithm_tables,
  path = "Results.xlsx"
)

# ============================================================
# 9. Record software environment for reproducibility
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)

