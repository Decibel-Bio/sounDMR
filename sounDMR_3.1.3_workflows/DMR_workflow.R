# Copyright 2023 Sound Agriculture Company
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# =============================================================================
# sounDMR workflow - CASE 1 (default): ONE target gene, shared by all treatments
# =============================================================================
#
# Use this when every treatment (oligo design, dose, ...) targets the same gene
# and you want to compare treatments against the control and each other.
#
# Set up for db210: T2-T5 (four oligo designs) all target EAB_chr6B and are each
# scored against the Water control. T6-T9 were designed for a different gene
# (EAB_chr6) that isn't scored in this experiment, so here they are used only as
# an off-target control: each of T2-T5 is also compared against the mean of T6-T9.
#
# Other designs:
#   CASE 2  one treatment targets several genes  -> DMR_workflow_multi_gene.R
#   CASE 3  treatments target different genes    -> DMR_workflow_treatment_specific_genes.R
#   (see README.md in this folder for a side-by-side of what changes)
#
#   Part 1  Load methyl beds, gene coordinates and experimental design   [same in all cases]
#   Part 2  Group DMR model                                             [same in all cases]
#   Part 3a Per treatment: changepoint regions + DMR scores             [same in all cases]
#   Part 3b Per treatment: bootstrap at THE target gene                 [CASE-SPECIFIC]
#   Part 4  Combine significant DMRs -> tables, BED, figures            [same call; read _summary.csv]
# =============================================================================

library(sounDMR)
set.seed(1)   # boot_score() samples random background positions


# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------

experiment_id <- "db210"                 # prefix for every output file
out_dir       <- paste0(experiment_id, "_downstream")

# Input files (use file.choose() instead of a path if you prefer to pick them)
bed_dir      <- "."                      # folder with the methyl bed files (18 for db210)
bed_pattern  <- "\\.methylbed$"          # e.g. 2026_10_09_155801_barcode01_merged.sorted_subset.methylbed
geneco_file  <- file.choose()            # gene coordinates: the target gene AND background genes
                                         #   (db210: Chromosome, Gene_Name, Strand, Low, High,
                                         #   Gene_ID, Adapt_Low, Adapt_High)
gene_column  <- "Gene_Name"              # gene-name column in the coordinates file
design_file  <- file.choose()            # design: ID, Bedfile, Plant, Group (C/T),
                                         #   Treatment (db210: Water, T2-T9)

# Experimental design columns / levels
target_variable <- "Treatment"           # fixed effect in the DMR model
control_label   <- "Water"               # control level of `target_variable`
control_group   <- "C"                   # control level of the Group column
treated_group   <- "T"                   # treated level of the Group column

# >>> CASE 1: the single target gene and the oligo window (relative to the gene, bp)
# db210: with dmr_score3 and a cutoff of 6 this window reproduces the saved
# DMR_boot_T2_sig ... DMR_boot_T5_sig tables (4 / 10 / 16 / 21 DMRs). With the
# defaults below (dmr_score4, final cutoff 8 in step 4.0) db210 gives 4 / 8 / 14 / 18.
target_gene  <- "EAB_chr6B"
target_start <- -150
target_end   <- 3000

# Treatments to score at target_gene (NULL = every non-control level).
# db210: only T2-T5 targeted EAB_chr6B; T6-T9 are the off-target control below.
treatments_to_score <- c("T2", "T3", "T4", "T5")

# Changepoint penalties (higher = fewer, larger regions)
cpt_penalty    <- c(CG = 3, CHG = 4, CHH = 6)   # db210: these reproduce the saved regions
save_cpt_plots <- TRUE

# Scoring
scoring_col    <- "dmr_score4"           # region score used by boot_score(): dmr_score, dmr_score2
                                         # (~10x smaller - lower min_sig_score), dmr_score3,
                                         # dmr_score4 (same scale as score3, but favours larger regions:
                                         #   x0.9 at n = 1, x1 at n = 5, x1.1 at n = 20, x1.2 at n = 100;
                                         #   0 if |mean Z| < 0.25 or |change| < 2),
                                         # dmr_score_noZ
min_sig_score  <- 6                      # adjusted_soundscore threshold for "significant" in Part 3
                                         #   (step 4.0 applies the final cutoff, 8 by default)

# Off-target treatments: treatments in this experiment that targeted a DIFFERENT
# gene (or a non-targeting / scrambled oligo). Each scored treatment is then also
# compared against their mean - <TX>_off per cytosine, Dif_Vs_Off_Target per region,
# and in <id>_dmrs.csv Percent_Dif_vs_Off_Target, Off_Target_Fraction and
# Off_Target_Adjusted_Score - alongside the control. NULL = control only.
# db210: T6-T9 targeted EAB_chr6, so they are the off-target control for T2-T5.
# (To also score off-target treatments at their own gene, use
#  DMR_workflow_treatment_specific_genes.R.)
off_target_trts <- c("T6", "T7", "T8", "T9")


# -----------------------------------------------------------------------------
# PART 1: Load data                                          [same in all cases]
# -----------------------------------------------------------------------------

methyl_bed <- list.files(path = bed_dir, pattern = bed_pattern, full.names = TRUE)
Geneco     <- read.csv(geneco_file, header = TRUE)

# With gene_info = TRUE this returns the zoomframe: only cytosines inside the
# windows in Geneco, with Zeroth_pos = position relative to the gene
Methylframe <- generate_methylframe(methyl_bed_list      = methyl_bed,
                                    Sample_count         = 0,
                                    Methyl_call_type     = "Modkit",
                                    filter_NAs           = 2,
                                    gene_info            = TRUE,
                                    gene_coordinate_file = Geneco,
                                    Gene_column          = gene_column,
                                    target_info          = FALSE,
                                    File_prefix          = "Sample")

# Needs ID (matching the bed sample IDs), Plant/Individual, Group and the
# `target_variable` column (e.g. Treatment = Water / T2 / T3 ...)
experimental_design_df <- read.csv(design_file, header = TRUE)
head(experimental_design_df)

if (is.null(treatments_to_score)) {
  treatments_to_score <- setdiff(unique(experimental_design_df[[target_variable]]), control_label)
}
stopifnot(control_label %in% experimental_design_df[[target_variable]],
          all(treatments_to_score %in% experimental_design_df[[target_variable]]),
          all(off_target_trts %in% experimental_design_df[[target_variable]]),
          target_gene %in% Geneco[[gene_column]])

# Every scored treatment uses the same off-target treatments (all hit target_gene)
off_targets <- setNames(rep(list(off_target_trts), length(treatments_to_score)),
                        treatments_to_score)


# -----------------------------------------------------------------------------
# PART 2: Group DMR model                                    [same in all cases]
# -----------------------------------------------------------------------------

dmr_obj <- create_dmr_obj(Methylframe, experimental_design_df)

# Per-cytosine summary. additional_summary_cols adds a <level>_mean column for
# every treatment level (T2_mean, Water_mean, ...), used below.
methyl_summary <- create_methyl_summary(dmr_obj,
                                        control = control_group,
                                        treated = treated_group,
                                        additional_summary_cols = list(c("mean", target_variable)))

# Binomial GLMM per cytosine: cbind(Meth, UnMeth) ~ (1|Individual) + Treatment,
# with the control as reference level. Adds Z_/Estimate_/StdErr_<Treatment><level>_small
# columns for every non-control level.
methyl_summary <- find_DMR(methyl_summary, dmr_obj,
                           fixed           = target_variable,
                           random          = "Individual",
                           reads_threshold = 3,
                           control         = control_label,
                           model           = "binomial",
                           analysis_type   = "group",
                           target_variable = target_variable)

all_genes <- unique(dmr_obj$ZoomFrame_filtered$Gene)


# -----------------------------------------------------------------------------
# PART 3a: Changepoint regions + DMR scores, once per treatment  [same in all cases]
# -----------------------------------------------------------------------------
# The model gives each treatment its own Z column, so regions and scores are
# always built per treatment, across all genes.

dmr_scores <- list()

for (trt in treatments_to_score) {
  z_col      <- paste0("Z_", target_variable, trt, "_small")
  mean_col   <- paste0(trt, "_mean")
  change_col <- paste0(trt, "_mean_dif")
  off_col    <- paste0(trt, "_off")
  if (!z_col %in% names(methyl_summary)) stop("No model column ", z_col, " - check the treatment name.")

  # treatment mean - control mean, per cytosine
  methyl_summary[[change_col]] <- methyl_summary[[mean_col]] - methyl_summary$Control

  # optional: treatment mean - mean of its off-target treatments (set in SETTINGS),
  # i.e. treatments that targeted a different gene
  off_trts <- off_targets[[trt]]
  if (length(off_trts)) {
    methyl_summary[[off_col]] <- methyl_summary[[mean_col]] -
      rowMeans(methyl_summary[, paste0(off_trts, "_mean"), drop = FALSE], na.rm = TRUE)
  }

  # split each gene x context into regions of similar Z score
  methyl_summary <- changepoint_analysis(methyl_summary,
                                         CG_penalty   = cpt_penalty[["CG"]],
                                         CHG_penalty  = cpt_penalty[["CHG"]],
                                         CHH_penalty  = cpt_penalty[["CHH"]],
                                         target_genes = all_genes,
                                         save_plots   = save_cpt_plots,
                                         z_col        = z_col)

  # summarise each region and give it a DMR score
  dmr_scores[[trt]] <- sound_score(changepoint_OF = methyl_summary,
                                   Statistic      = z_col,
                                   Per_Change     = change_col,
                                   other_columns  = if (length(off_trts)) off_col else character(0),
                                   CF = FALSE, UserFunction = NA)
}


# -----------------------------------------------------------------------------
# PART 3b: Bootstrap at the target gene, once per treatment     [CASE 1]
# -----------------------------------------------------------------------------
# Every treatment is scored at the same gene and window. All other genes are the
# bootstrap background. One sig table per treatment: DMR_boot_<TX>_sig.

dmr_boots <- list()
sig_list  <- list()

for (trt in treatments_to_score) {
  key <- paste0("DMR_boot_", trt)
  dmr_boots[[key]] <- boot_score(sound_score_obj  = dmr_scores[[trt]],
                                 target_gene      = target_gene,
                                 scoring_col_name = scoring_col,
                                 target_start     = target_start,
                                 target_end       = target_end)
  sig_list[[paste0(key, "_sig")]] <- get_sig_dmrs(dmr_boots[[key]], min_score = min_sig_score)
}


# -----------------------------------------------------------------------------
# PART 4: Combine significant DMRs -> tables, BED and figures
# -----------------------------------------------------------------------------
# CASE 1: <id>_summary.csv (one row per treatment) is the main table;
# <id>_summary_by_gene.csv is identical apart from the Gene column.

# -----------------------------------------------------------------------------
# 4.0  Final significance cutoff, applied before building the tables
#      (no need to rerun Part 3 - the bootstrap results are all still in dmr_boots)
# -----------------------------------------------------------------------------
# One cutoff for every treatment:        new_min_sig_score <- 8   (default)
# Or per treatment (others unchanged):   new_min_sig_score <- c(T2 = 7, T4 = 5)
# NULL = keep min_sig_score from SETTINGS
new_min_sig_score <- 8

if (!is.null(new_min_sig_score)) {
  sig_trts <- treatments_to_score   # treatment for each element of sig_list
  for (i in seq_along(sig_list)) {
    t <- sig_trts[i]
    if (!is.null(names(new_min_sig_score)) && !t %in% names(new_min_sig_score)) next
    cutoff <- if (is.null(names(new_min_sig_score))) new_min_sig_score else new_min_sig_score[[t]]
    sig_list[[i]] <- get_sig_dmrs(dmr_boots[[sub("_sig$", "", names(sig_list)[i])]], min_score = cutoff)
    cat(names(sig_list)[i], ": cutoff", cutoff, "->", nrow(sig_list[[i]]), "DMRs\n")
  }
}

res <- run_dmr_downstream(sig_list,
                          dmr_obj,
                          control         = control_label,
                          treatments      = treatments_to_score,   # order matches sig_list
                          off_targets     = off_targets,       # off-target % methylation per DMR
                          target_variable = target_variable,
                          out_dir         = out_dir,
                          prefix          = experiment_id)
                          # design = read.csv("oligo_info.csv")  # optional: fills Description /
                          #   Oligo_Duplexes (columns: Treatment + annotation)

res$summary

# Combined methyl_summary (every treatment's columns) -> <id>_ms.csv
data.table::fwrite(methyl_summary, file.path(out_dir, paste0(experiment_id, "_ms.csv")))

# Whole session -> <id>.RData (methyl_summary, dmr_obj, dmr_scores, dmr_boots, sig_list,
# res, ...). Reload later with load("<out_dir>/<id>.RData").
save.image(file.path(out_dir, paste0(experiment_id, ".RData")))
