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
# sounDMR workflow - STEPWISE: one treatment at a time (works for cases 1, 2, 3)
#   >>> filled in for db210: T2-T5 target EAB_chr6B, T6-T9 (EAB_chr6) are their
#       off-targets, Water control
#       Run Parts 1-2 once. Run Part 3 with trt <- "T2", then change it to "T3" and run
#       Part 3 again, then "T4", then "T5". Then run Part 4.
# =============================================================================
#
# Same steps as the looped scripts, but Part 3 is a block you run by hand for one
# treatment at a time. For each treatment you set its changepoint penalties, target
# gene(s) + window, score direction and significance cutoff, check the results,
# rerun with different settings if needed, then move on to the next treatment.
#
#   Part 1  Load methyl beds, gene coordinates and experimental design   [run once]
#   Part 2  Group DMR model                                             [run once]
#   Part 3  ONE TREATMENT: 3.1 settings -> 3.2 changepoints -> 3.3 region scores
#           -> 3.4 bootstrap -> 3.5 cutoff + check -> 3.6 keep            [repeat per treatment]
#   Part 4  Combine everything kept -> tables, BED, figures, settings log [run once at the end]
#
# How the three designs fit in (step 3.1, `targets_trt`):
#   CASE 1  same gene for every treatment   -> one row, the same gene in every block
#   CASE 2  a treatment targets several genes -> one row per gene in that block
#   CASE 3  treatments target different genes -> one row, a different gene per block
#
# All treatments share one methyl_summary: each block adds that treatment's
# columns (<TX>_mean_dif, MethRegion_Z_<target_variable><TX>_small, ...) to it.
# Rerunning a block for the same treatment overwrites its columns and replaces
# what was kept for it, so you can iterate freely.
# =============================================================================

library(sounDMR)
set.seed(1)   # boot_score() samples random background positions


# -----------------------------------------------------------------------------
# SETTINGS (shared by all treatments; per-treatment settings are in step 3.1)
# -----------------------------------------------------------------------------

experiment_id <- "db210"                 # prefix for every output file
out_dir       <- paste0(experiment_id, "_downstream")

# Input files (use file.choose() instead of a path if you prefer to pick them)
bed_dir      <- "."                      # folder with the 18 db210 methyl beds
bed_pattern  <- "\\.methylbed$"          # e.g. 2026_10_09_155801_barcode01_merged.sorted_subset.methylbed
geneco_file  <- file.choose()            # db210 gene coordinates (Chromosome, Gene_Name, Strand,
                                         #   Low, High, Gene_ID, Adapt_Low, Adapt_High)
gene_column  <- "Gene_Name"
design_file  <- file.choose()            # db210 design: ID (S1-S18), Bedfile, Plant, Group (C/T),
                                         #   Treatment (Water, T2-T9)

# Experimental design columns / levels
target_variable <- "Treatment"           # fixed effect in the DMR model
control_label   <- "Water"               # control level of `target_variable`
control_group   <- "C"                   # control level of the Group column
treated_group   <- "T"                   # treated level of the Group column

# Region score that boot_score() distance-adjusts into adjusted_soundscore (and that
# min_sig_score is applied to). Kept the same for all treatments so their summaries
# are comparable. Options from sound_score():
#   "dmr_score"     n^(1/3) * |mean Z| * |mean % change|^(1/2)
#   "dmr_score2"    like dmr_score, but uses the arcsine-transformed change
#                   -> roughly 10x smaller values; lower min_sig_score to match
#   "dmr_score3"    n^(1/3) * |mean Z|^(1/2) * |mean % change|^(1/2) * |max Z|^(1/2)
#   "dmr_score4"    dmr_score3 with n^(1/2.5) instead of n^(1/3): slightly favours larger
#                   regions (x1.17 at n = 10, x1.36 at n = 100), but only for regions with a
#                   real average shift: 0 if |mean Z| < 0.25 or |mean % change| < 2, so long
#                   flat regions can't score on size alone. Tune in sound_score():
#                   score4_count_exponent, score4_min_mean_z, score4_min_mean_dif
#   "dmr_score_noZ" n^(1/3) * |mean % change|^(1/2)   (ignores the model)
scoring_col <- "dmr_score3"


# -----------------------------------------------------------------------------
# PART 1: Load data                                          [run once]  
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
stopifnot(control_label %in% experimental_design_df[[target_variable]])


# -----------------------------------------------------------------------------
# PART 2: Group DMR model                                    [run once]  
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

# Results collected across treatment blocks
dmr_scores     <- list()        # sound_score() object per treatment
dmr_boots      <- list()        # boot_score() object per treatment x gene
sig_list       <- list()        # kept significant-DMR tables
sig_treatments <- character(0)  # treatment for each element of sig_list (same names)
settings_log   <- list()        # settings used for each treatment x gene
off_targets    <- list()        # off-target treatments chosen for each treatment


# =============================================================================
# PART 3: ONE TREATMENT AT A TIME
#   Run 3.1 -> 3.6, check, adjust 3.1 and rerun if needed. Then set the next
#   treatment in 3.1 and run 3.1 -> 3.6 again.
# =============================================================================

# -----------------------------------------------------------------------------
# 3.1  Choose the treatment and its settings
# -----------------------------------------------------------------------------
trt <- "T2"                    # <<< change to "T3", "T4", "T5" and rerun 3.1 -> 3.6

# Target gene(s) for THIS treatment, with the oligo window relative to each gene (bp).
#   one row     -> cases 1 and 3
#   several rows -> case 2 (the other rows are left out of each gene's bootstrap background)
# db210: every treatment targets EAB_chr6B. This window, with dmr_score3 and a cutoff
# of 6, reproduces the saved DMR_boot_T2_sig ... DMR_boot_T5_sig tables.
targets_trt <- data.frame(Gene = "EAB_chr6B", target_start = -150, target_end = 3000)

cpt_penalty    <- c(CG = 3, CHG = 4, CHH = 6)  # higher = fewer, larger regions
                                               # (these reproduce the saved db210 regions)
direction      <- "positive"                   # "positive" (gain of methylation), "negative", "absolute"
min_sig_score  <- 6                            # adjusted_soundscore cutoff

# Off-target treatments: T6-T9 targeted EAB_chr6, not EAB_chr6B, so they are the
# off-target baseline for T2-T5 (same for every block). This adds <TX>_off and, in
# _dmrs.csv, Off_Target_Fraction and Off_Target_Adjusted_Score (the last two columns
# of the db210 DMRs sheet). Set to NULL to skip.
off_target_trts <- c("T6", "T7", "T8", "T9")

z_col      <- paste0("Z_", target_variable, trt, "_small")
mean_col   <- paste0(trt, "_mean")
change_col <- paste0(trt, "_mean_dif")
off_col    <- paste0(trt, "_off")
stopifnot(z_col %in% names(methyl_summary),
          all(targets_trt$Gene %in% methyl_summary$Gene),
          !trt %in% off_target_trts,
          length(off_target_trts) == 0 ||
            all(paste0(off_target_trts, "_mean") %in% names(methyl_summary)))

# -----------------------------------------------------------------------------
# 3.2  Changepoint regions for this treatment
# -----------------------------------------------------------------------------
# treatment mean - control mean, per cytosine
methyl_summary[[change_col]] <- methyl_summary[[mean_col]] - methyl_summary$Control

# optional: treatment mean - mean of its off-target treatments, per cytosine
if (length(off_target_trts)) {
  methyl_summary[[off_col]] <- methyl_summary[[mean_col]] -
    rowMeans(methyl_summary[, paste0(off_target_trts, "_mean"), drop = FALSE], na.rm = TRUE)
}

# Regions are found for every gene; plots are only saved for this block's target genes.
# Adds this treatment's MethRegion_/MethGroup_ columns to the shared methyl_summary.
methyl_summary <- changepoint_analysis(methyl_summary,
                                       CG_penalty   = cpt_penalty[["CG"]],
                                       CHG_penalty  = cpt_penalty[["CHG"]],
                                       CHH_penalty  = cpt_penalty[["CHH"]],
                                       target_genes = targets_trt$Gene,
                                       save_plots   = TRUE,
                                       z_col        = z_col)

# CHECK: open <gene>_<CX>_changepoints_line.png / _volcano.png in the working directory.
#   Regions too fragmented -> raise that context's penalty in 3.1 and rerun 3.2.
#   Real signal merged into large flat regions -> lower it.
#   (The PNGs are overwritten by the next treatment, so check them now.)

# -----------------------------------------------------------------------------
# 3.3  Score every region
# -----------------------------------------------------------------------------
dmr_scores[[trt]] <- sound_score(changepoint_OF = methyl_summary,
                                 Statistic      = z_col,
                                 Per_Change     = change_col,
                                 other_columns  = if (length(off_target_trts)) off_col else character(0),
                                 CF = FALSE, UserFunction = NA)

# CHECK: strongest regions on the target gene(s)
rs_trt <- dmr_scores[[trt]]$region_summary
rs_trt <- rs_trt[rs_trt$Gene %in% targets_trt$Gene, ]
head(rs_trt[order(-rs_trt[[scoring_col]]),
            c("cp_group", "CX", "Start", "Stop", "Count", paste0(change_col, "_mean"), scoring_col)], 15)

# -----------------------------------------------------------------------------
# 3.4  Bootstrap at each target gene
# -----------------------------------------------------------------------------
boot_keys <- paste0("DMR_boot_", trt, "_", targets_trt$Gene)
for (i in seq_len(nrow(targets_trt))) {
  gene <- targets_trt$Gene[i]
  dmr_boots[[boot_keys[i]]] <- boot_score(sound_score_obj  = dmr_scores[[trt]],
                                          target_gene      = gene,
                                          scoring_col_name = scoring_col,
                                          target_start     = targets_trt$target_start[i],
                                          target_end       = targets_trt$target_end[i],
                                          direction_DMR    = direction,
                                          exclude_genes    = setdiff(targets_trt$Gene, gene))
}

# CHECK: the printed bootstrap p-values and "Final Bootstrap Adjusted ... Score" per
# context (<1 nothing, 1-2 subtle, 2-3 moderate, 3-5 significant, 5+ very strong).

# -----------------------------------------------------------------------------
# 3.5  Apply the cutoff and look at what passes
# -----------------------------------------------------------------------------
sig_trt <- lapply(dmr_boots[boot_keys], get_sig_dmrs, min_score = min_sig_score)
names(sig_trt) <- paste0(boot_keys, "_sig")

# CHECK: how many pass, by gene and context, and how the scores are distributed
sapply(sig_trt, nrow)
lapply(sig_trt, function(s) table(s$CX))
lapply(dmr_boots[boot_keys], function(b) summary(b$target_rs$adjusted_soundscore))

# CHECK: where the passing regions sit relative to the target window
if (sum(sapply(sig_trt, nrow)) > 0) {
  print(plot_dmr_map(combine_sig_dmrs(sig_trt, treatments = rep(trt, length(sig_trt)),
                                      verbose = FALSE)))
}
# Too many / too few -> change min_sig_score in 3.1 and rerun 3.5 (no need to redo 3.2-3.4).

# -----------------------------------------------------------------------------
# 3.6  Keep this treatment's results (replaces anything kept earlier for it)
# -----------------------------------------------------------------------------
keep           <- sig_treatments != trt
sig_list       <- c(sig_list[keep], sig_trt)
sig_treatments <- c(sig_treatments[keep], setNames(rep(trt, length(sig_trt)), names(sig_trt)))
off_targets[[trt]] <- off_target_trts           # NULL removes any earlier entry
settings_log[[trt]] <- data.frame(
  Treatment = trt, Gene = targets_trt$Gene,
  target_start = targets_trt$target_start, target_end = targets_trt$target_end,
  CG_penalty = cpt_penalty[["CG"]], CHG_penalty = cpt_penalty[["CHG"]],
  CHH_penalty = cpt_penalty[["CHH"]], direction = direction,
  min_sig_score = min_sig_score,
  off_targets = if (length(off_target_trts)) paste(off_target_trts, collapse = "; ") else NA,
  n_sig = sapply(sig_trt, nrow), row.names = NULL)

settings_log[[trt]]

# >>> Next treatment: change `trt` (and any settings) in 3.1 and run 3.1 -> 3.6 again.


# =============================================================================
# PART 4: Combine every kept treatment -> tables, BED and figures
# =============================================================================
# Main table: <id>_summary.csv if every treatment has one gene (case 1),
#             <id>_summary_by_gene.csv if any treatment has several genes (cases 2/3).

# -----------------------------------------------------------------------------
# 4.0  OPTIONAL: change the significance cutoff for kept treatments (no need to
#      rerun Part 3 - the bootstrap results are all still in dmr_boots)
# -----------------------------------------------------------------------------
# One cutoff for every treatment:        new_min_sig_score <- 7
# Or per treatment (others unchanged):   new_min_sig_score <- c(T2 = 7, T4 = 5)
new_min_sig_score <- NULL      # NULL = keep the cutoffs chosen in each block

if (!is.null(new_min_sig_score)) {
  for (k in names(sig_list)) {
    t <- sig_treatments[[k]]
    if (!is.null(names(new_min_sig_score)) && !t %in% names(new_min_sig_score)) next
    cutoff <- if (is.null(names(new_min_sig_score))) new_min_sig_score else new_min_sig_score[[t]]
    sig_list[[k]] <- get_sig_dmrs(dmr_boots[[sub("_sig$", "", k)]], min_score = cutoff)
    # keep the settings log in step with what was actually used
    row <- paste0("DMR_boot_", t, "_", settings_log[[t]]$Gene, "_sig") == k
    settings_log[[t]]$min_sig_score[row] <- cutoff
    settings_log[[t]]$n_sig[row]         <- nrow(sig_list[[k]])
  }
  print(do.call(rbind, settings_log)[, c("Treatment", "Gene", "min_sig_score", "n_sig")])
}

res <- run_dmr_downstream(sig_list,
                          dmr_obj,
                          control         = control_label,
                          treatments      = unname(sig_treatments),   # order matches sig_list
                          off_targets     = off_targets,              # per-treatment off-target sets
                          target_variable = target_variable,
                          out_dir         = out_dir,
                          prefix          = experiment_id)
                          # design = read.csv("oligo_info.csv")  # optional: fills Description /
                          #   Oligo_Duplexes (columns: Treatment [+ Gene] + annotation)

res$summary
res$summary_by_gene

# Record the treatment-specific settings next to the results
settings_log_df <- do.call(rbind, settings_log)
write.csv(settings_log_df, file.path(out_dir, paste0(experiment_id, "_settings_log.csv")),
          row.names = FALSE)

# Combined methyl_summary (every treatment's columns) -> <id>_ms.csv
data.table::fwrite(methyl_summary, file.path(out_dir, paste0(experiment_id, "_ms.csv")))

# Whole session -> <id>.RData (methyl_summary, dmr_obj, dmr_scores, dmr_boots, sig_list,
# res, ...). Reload later with load("<out_dir>/<id>.RData").
save.image(file.path(out_dir, paste0(experiment_id, ".RData")))
