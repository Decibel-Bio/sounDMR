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
# sounDMR workflow - CASE 2: ONE treatment targets SEVERAL genes
# =============================================================================
#
# Use this when a treatment (e.g. a pooled oligo set) was designed to hit more
# than one gene, and every listed treatment shares the same set of target genes.
#   e.g. db187: T2_PAL2 targets LsPAL1 and LsPAL2
#
# Other designs:
#   CASE 1  one target gene for all treatments  -> DMR_workflow.R (default)
#   CASE 3  treatments target different genes    -> DMR_workflow_treatment_specific_genes.R
#   (see README.md in this folder for a side-by-side of what changes)
#
#   Part 1  Load methyl beds, gene coordinates and experimental design   [same in all cases]
#   Part 2  Group DMR model                                             [same in all cases]
#   Part 3a Per treatment: changepoint regions + DMR scores             [same in all cases]
#           -> computed ONCE per treatment and reused for all its genes
#   Part 3b Per treatment x gene: bootstrap, with the treatment's OTHER
#           target genes left out of the background                     [CASE-SPECIFIC]
#   Part 4  Combine significant DMRs -> tables, BED, figures            [same call; read _summary_by_gene.csv]
# =============================================================================

library(sounDMR)
set.seed(1)   # boot_score() samples random background positions


# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------

experiment_id <- "db187"                 # prefix for every output file
out_dir       <- paste0(experiment_id, "_downstream")

# Input files (use file.choose() instead of a path if you prefer to pick them)
bed_dir      <- "."                      # folder with the methyl bed files
bed_pattern  <- "\\.bed$"
geneco_file  <- "Geneco.csv"             # gene coordinates: the target gene AND background genes
gene_column  <- "Gene_Name"              # gene-name column in the coordinates file
design_file  <- "experimental_design.csv"

# Experimental design columns / levels
target_variable <- "Treatment"           # fixed effect in the DMR model
control_label   <- "Water"               # control level of `target_variable`
control_group   <- "C"                   # control level of the Group column
treated_group   <- "T"                   # treated level of the Group column

# >>> CASE 2: the treatment(s), and the genes each of them targets (one row per gene)
#   Gene              - target gene name, as it appears in the Gene column
#   target_start/_end - oligo window relative to that gene (bp)
treatments_to_score <- c("T2_PAL2")      # every treatment listed here is scored at every gene below
target_genes <- data.frame(
  Gene         = c("LsPAL1", "LsPAL2"),
  target_start = c(-2000,    -2000),
  target_end   = c( 1000,     1000)
)

# Leave the treatment's other target genes out of each gene's bootstrap
# background (their real DMRs would otherwise inflate it -> conservative scores)
exclude_other_targets <- TRUE

# Changepoint penalties (higher = fewer, larger regions)
cpt_penalty    <- c(CG = 3, CHG = 4, CHH = 6)
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

# Off-target treatments (optional). Treatments in this experiment that targeted a
# DIFFERENT gene (or a non-targeting / scrambled oligo). Each scored treatment is
# then also compared against their mean - <TX>_off per cytosine, Dif_Vs_Off_Target
# per region, Percent_Dif_vs_Off_Target per DMR - alongside the control.
#   e.g. T2_PAL2 targeted LsPAL1 + LsPAL2 and T5, T6 targeted other genes:
#   off_target_trts <- c("T5", "T6")
# (To also score T5/T6 at their own genes, use DMR_workflow_treatment_specific_genes.R.)
off_target_trts <- NULL


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

stopifnot(control_label %in% experimental_design_df[[target_variable]],
          all(treatments_to_score %in% experimental_design_df[[target_variable]]),
          all(off_target_trts %in% experimental_design_df[[target_variable]]),
          all(target_genes$Gene %in% Geneco[[gene_column]]))

# Every scored treatment uses the same off-target treatments
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
# PART 3b: Bootstrap per treatment x target gene                [CASE 2]
# -----------------------------------------------------------------------------
# The same sound_score object (one per treatment) is reused for every gene; only
# the target gene / window changes. With exclude_other_targets = TRUE the
# treatment's other target genes are not used as background.
# One sig table per treatment x gene: DMR_boot_<TX>_<gene>_sig.

dmr_boots      <- list()
sig_list       <- list()
sig_treatments <- character(0)

for (trt in treatments_to_score) {
  for (i in seq_len(nrow(target_genes))) {
    gene <- target_genes$Gene[i]
    key  <- paste0("DMR_boot_", trt, "_", gene)
    dmr_boots[[key]] <- boot_score(sound_score_obj  = dmr_scores[[trt]],
                                   target_gene      = gene,
                                   scoring_col_name = scoring_col,
                                   target_start     = target_genes$target_start[i],
                                   target_end       = target_genes$target_end[i],
                                   exclude_genes    = if (exclude_other_targets)
                                                        setdiff(target_genes$Gene, gene))
    sig_list[[paste0(key, "_sig")]] <- get_sig_dmrs(dmr_boots[[key]], min_score = min_sig_score)
    sig_treatments <- c(sig_treatments, trt)
  }
}


# -----------------------------------------------------------------------------
# PART 4: Combine significant DMRs -> tables, BED and figures
# -----------------------------------------------------------------------------
# CASES 2/3: <id>_summary_by_gene.csv (one row per treatment x gene) is the main
# table. <id>_summary.csv pools all genes per treatment and lists them in `Genes`.
# The DMR map PDF has one column of panels per gene.

# -----------------------------------------------------------------------------
# 4.0  Final significance cutoff, applied before building the tables
#      (no need to rerun Part 3 - the bootstrap results are all still in dmr_boots)
# -----------------------------------------------------------------------------
# One cutoff for every treatment:        new_min_sig_score <- 8   (default)
# Or per treatment (others unchanged):   new_min_sig_score <- c(T2 = 7, T4 = 5)
# NULL = keep min_sig_score from SETTINGS
new_min_sig_score <- 8

if (!is.null(new_min_sig_score)) {
  sig_trts <- sig_treatments   # treatment for each element of sig_list
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
                          treatments      = sig_treatments,        # order matches sig_list
                          off_targets     = off_targets,       # off-target % methylation per DMR
                          target_variable = target_variable,
                          out_dir         = out_dir,
                          prefix          = experiment_id)
                          # design = read.csv("oligo_info.csv")  # optional: fills Description /
                          #   Oligo_Duplexes (columns: Treatment, Gene + annotation)

res$summary_by_gene

# Combined methyl_summary (every treatment's columns) -> <id>_ms.csv
data.table::fwrite(methyl_summary, file.path(out_dir, paste0(experiment_id, "_ms.csv")))

# Whole session -> <id>.RData (methyl_summary, dmr_obj, dmr_scores, dmr_boots, sig_list,
# res, ...). Reload later with load("<out_dir>/<id>.RData").
save.image(file.path(out_dir, paste0(experiment_id, ".RData")))
