#!/usr/bin/env Rscript
#
# Purpose:
# R script to summarize sequencing read counts obtained from a CRISPRi library.
# The script summarizes count tables per sample into one main table, adds
# statistical metrics for a pairwise sample comparison using DESeq2, and calculates
# fitness scores for each gene and condition.
#
# Date: 2024-07-05
# Author: Ute Hoffmann
# Affilation: Science For Life Laboratory (KTH), Stockholm, Sweden
# based on work by Michael Jahn, PhD

# PARSE INPUT
# ====================
#
# input arguments
setwd("~/Documents/bioinformatic_analyses/2024-05_LargeRubiscoLibrary/nonWeighted_fitnessCalc")
path_samplesheet <- "../input/samplesheet_1stCultiv_reduced.csv" # file paths to sample sheet
path_counts <- "../output_selectedBarcodes/prepare/all_counts.tsv" # file paths to count tables
normalization <- FALSE # default: FALSE
gene_fitness <- TRUE # default: TRUE
gene_sep <- "|" # default: "|" aka the pipe symbol
gene_controls <- "" # default: "" aka empty string
number_cores <- 4 # number of CPU cores


# LOAD PACKAGES
# ====================
#
library(readr)
library(tibble)
library(stringr)
library(tidyr)
library(dplyr)
library(purrr)
library(DESeq2)
library(BiocParallel)

if (normalization) {
    library(limma)
}

# DATA PREPARATION
# ====================
#
# Step 1: Load sample layout sheet - file names must be row names
df_samplesheet <- readr::read_csv(path_samplesheet, col_types = cols()) %>%
    select(all_of(c("sample", "condition", "replicate", "time", "group", "reference_group"))) %>%
    dplyr::mutate(group = factor(`group`))

# Check which conditions
# - do not have differing groups/reference when zero time
# - have the minimum of two distinct time points OR
# - are compared to the zero time point of another condition OR
# - for all others skip DESeq2 and fitness calculation
df_samplesheet <- df_samplesheet %>%
    dplyr::group_by(condition) %>%
    filter(!(time == 0 & group != reference_group)) %>%
    filter(length(unique(time)) >= 2 || (time != 0 & group != reference_group)) %>%
    tibble::column_to_rownames("sample")

df_counts <- readr::read_tsv(path_counts, col_types = cols())
df_counts <- tidyr::pivot_longer(df_counts,
    cols = 3:ncol(df_counts),
    names_to = "sample", values_to = "numreads"
)

# print overview information to console
message("Number of sgRNAs detected in n samples:")
df_counts %>%
    group_by(sgRNA) %>%
    dplyr::summarize(sgRNAs_detected_in_samples = sum(numreads > 0)) %>%
    dplyr::count(sgRNAs_detected_in_samples) %>%
    dplyr::mutate(percent_total = n / sum(n) * 100) %>%
    dplyr::arrange(dplyr::desc(sgRNAs_detected_in_samples)) %>%
    print(n=80)

# input data frame must be reshaped to a 'counts matrix' with genes as rows
# and samples (conditions) as columns.
counts <- df_counts %>%
    dplyr::select(-Gene) %>%
    # spread condition over columns and sgRNAs over rows
    tidyr::pivot_wider(names_from = sample, values_from = numreads) %>%
    # remove sgRNA column, replace NA with 0
    dplyr::mutate_at(vars(-1), function(x) coalesce(x, 0)) %>%
    # add row_names from column
    tibble::column_to_rownames("sgRNA")

# DIFFERENTIAL ABUNDANCE
# ======================
#
# DESeq2 can be used to obtain fold changes and significance metrics
# for condition-wise comparisons, for details see publication:
# Love, M.I., Huber, W., Anders, S. Genome Biology, 15:550, 2014.
# (https://doi.org/10.1186/s13059-014-0550-8)

message("Running DESeq2 for pairwise comparison.\nNote: this step can be time and computation-intense.")
message(paste0("Number of CPU cores used for DESeq parallel execution: ", number_cores, "."))

# 2. Meta data
# Meta data is required to carry out the actual DESeq2 analysis
# by 'contrasting' (comparing) selected conditions to each other.
# We check that the order of file names corresponds to colnames of counts
if (!all(row.names(df_samplesheet) %in% colnames(counts))) {
    stop("Not all samples listed in the sample sheet have
        corresponding read count data.")
}
counts <- counts[row.names(df_samplesheet)]

# 3. Perform DESeq2 analysis
DESeq_result <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = df_samplesheet,
    design = ~group
) %>%
    DESeq()

# get normalized count data
norm_counts <- counts(DESeq_result, normalized=TRUE)
t_norm_counts <- as.data.frame(t(norm_counts))
t_norm_counts$samplename <- row.names(t_norm_counts)
df_samplesheet$samplename <- row.names(df_samplesheet)
annot_norm_counts <- left_join(t_norm_counts, df_samplesheet)
annot_norm_counts <- subset(annot_norm_counts, annot_norm_counts$replicate < 8)
annot_norm_counts_long <- annot_norm_counts %>% pivot_longer(!c("samplename", "condition", "replicate", "time", "group", "reference_group"), names_to="sgRNA", values_to="count")
annot_norm_counts_long$samplename <- NULL
annot_norm_counts_long$group <- NULL
annot_norm_counts_long$reference_group <- NULL
annot_norm_counts_wide <- annot_norm_counts_long %>% pivot_wider(names_from="time", values_from="count")
annot_norm_counts_wide

calc_auc <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
}
time_points <- sort(unique(annot_norm_counts_long$time))[c(2,3,4)]

log2_transform <- annot_norm_counts_wide %>% rowwise() %>% mutate(`3.7`=log2(`3.7`/`0`), `8.6`=log2(`8.6`/`0`), `10.1`=log2(`10.1`/`0`))
subs_log2_transform <- subset(log2_transform, is.finite(as.vector(log2_transform[,5])[[1]]))
subs_log2_transform <- subset(subs_log2_transform, is.finite(as.vector(subs_log2_transform[,6])[[1]]))
subs_log2_transform <- subset(subs_log2_transform, is.finite(as.vector(subs_log2_transform[,7])[[1]]))

a <- subs_log2_transform %>% rowwise() %>% mutate(auc=calc_auc(time_points, c(`3.7`, `8.6`, `10.1`))/max(time_points)/2)
a$`0` <- NULL
a$`3.7` <- NULL
a$`8.6` <- NULL
a$`10.1` <- NULL
a <- pivot_wider(a, names_from=replicate, values_from=auc)
a <- a %>% rowwise() %>% mutate(mean_auc=mean(c(`1`,`2`,`3`,`4`,`5`,`6`,`7`), na.rm=TRUE), sd_auc=sd(c(`1`,`2`,`3`,`4`,`5`,`6`,`7`), na.rm=TRUE))
a[,c(3:9)] <- NULL


b <- a %>%
  # split sgRNA names into target gene and position
  tidyr::separate(sgRNA,
                  into = c("sgRNA_target", "sgRNA_position"), sep = paste0("\\", gene_sep),
                  remove = FALSE
  )

b$weight_auc <- b$mean_auc/(b$sd_auc^2)

b[is.na(as.vector(b$sd_auc)),]$sd_auc <- max(b[!is.na(as.vector(b$sd_auc)),]$sd_auc)

DESeq_result_table <- b
DESeq_result_table <- dplyr::left_join(
  DESeq_result_table,
  DESeq_result_table %>%
    dplyr::group_by(sgRNA_target, condition) %>%
    dplyr::summarize(
      .groups = "keep",
      # gene fitness
      wmean_fitness = weighted.mean(mean_auc, 1/(sd_auc * sd_auc), na.rm=TRUE), 
      sd_fitness = sd(mean_auc)
    ),
  by = c("sgRNA_target", "condition")
)

WT_fitness <- unique(DESeq_result_table[DESeq_result_table$sgRNA_target=="WT",]$wmean_fitness)
WT_fitness
K214_fitness <- max(unique(subset(DESeq_result_table, grepl("K214", DESeq_result_table$sgRNA_target))$wmean_fitness))
K214_fitness
DESeq_result_table$norm_fit <- (DESeq_result_table$wmean_fitness-K214_fitness)/(WT_fitness-K214_fitness)

p <- ggplot(DESeq_result_table, aes(x=norm_fit)) + geom_histogram(bins=300)
p

# EXPORT PROCESSED DATA
# =====================
#
# Save result tables to output folder, in Rdata format
if (nrow(df_samplesheet) > 0) {
    message("Saving 'result_Weighted_AUC_sd.Rdata'")
    save(DESeq_result_table, file = "result_Weighted_AUC_sd.Rdata")
    message("Saving 'result_Weighted_AUC_sd.tsv'")
    readr::write_tsv(DESeq_result_table, file = "result_Weighted_AUC_sd.tsv")
}

