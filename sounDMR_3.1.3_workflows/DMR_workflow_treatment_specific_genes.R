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
# sounDMR workflow - CASE 3: DIFFERENT treatments target DIFFERENT genes
# =============================================================================
#
# Use this when each treatment has its own target gene(s).
#   e.g. T2, T3, T4 -> LsPAL1 and T5, T6, T7 -> LsPAL2
# A treatment may have several rows (several genes), so this also covers mixes
# of CASE 2 and CASE 3.
#
# Other designs:
#   CASE 1  one target gene for all treatments  -> DMR_workflow.R (default)
#   CASE 2  one treatment targets several genes  -> DMR_workflow_multi_gene.R
#   (see README.md in this folder for a side-by-side of what changes)
#
#   Part 1  Load methyl beds, gene coordinates and experimental design   [same in all cases]
#   Part 2  Group DMR model                                             [same in all cases]
#   Part 3a Per treatment: changepoint regions + DMR scores             [same in all cases]
#   Part 3b Per row of `targets` (treatment x gene): bootstrap at that
#           treatment's own gene                                        [CASE-SPECIFIC]
#   Part 4  Combine significant DMRs -> tables, BED, figures            [same call; read _summary_by_gene.csv]
# =============================================================================

library(sounDMR)
set.seed(1)   # boot_score() samples random background positions


# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------

experiment_id <- "dbXXX"                 # prefix for every output file
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

# >>> CASE 3: one row per treatment x target gene
#   Treatment         - level of `target_variable`
#   Gene              - that treatment's target gene, as it appears in the Gene column
#   target_start/_end - oligo window relative to that gene (bp)
targets <- data.frame(
  Treatment    = c("T2",     "T3",     "T4",     "T5",     "T6",     "T7"),
  Gene         = c("LsPAL1", "LsPAL1", "LsPAL1", "LsPAL2", "LsPAL2", "LsPAL2"),
  target_start = c(-2000,    -2000,    -2000,    -1500,    -1500,    -1500),
  target_end   = c( 1000,     1000,     1000,      500,      500,      500)
)
# A treatment can have several rows (several genes), e.g. add
#   T7 / LsAGL24 / -1000 / 0   if T7 also targeted LsAGL24

# If a treatment has several target genes, leave its OTHER target genes out of
# each gene's bootstrap background (no effect for one-gene treatments).
# Genes targeted only by other treatments stay in the background: this
# treatment's Z scores there should be null.
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

# Off-target treatments: for each treatment, the treatments that targeted a
# DIFFERENT gene. Each treatment is then also compared against their mean - <TX>_off
# per cytosine, Dif_Vs_Off_Target per region, Percent_Dif_vs_Off_Target per DMR -
# alongside the control.
# Built automatically from `targets`: every other treatment that shares none of this
# treatment's genes. With the example above:
#   T2, T3, T4 -> off-targets T5, T6, T7     T5, T6, T7 -> off-targets T2, T3, T4
off_targets <- sapply(unique(targets$Treatment), function(t) {
  shares_a_gene <- unique(targets$Treatment[targets$Gene %in% targets$Gene[targets$Treatment == t]])
  setdiff(unique(targets$Treatment), shares_a_gene)
}, simplify = FALSE)
# Edit by hand if needed, e.g. add a scrambled-oligo treatment that isn't in `targets`:
#   off_targets$T2 <- c("T5", "T6", "T7", "T9")
# or skip the off-target comparison entirely:
#   off_targets <- list()


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

treatments_to_score <- unique(targets$Treatment)
stopifnot(control_label %in% experimental_design_df[[target_variable]],
          all(treatments_to_score %in% experimental_design_df[[target_variable]]),
          all(unlist(off_targets) %in% experimental_design_df[[target_variable]]),
          all(targets$Gene %in% Geneco[[gene_column]]))

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
# PART 3b: Bootstrap per row of `targets`                      [CASE 3]
# -----------------------------------------------------------------------------
# Each treatment is scored only at its own gene(s), using that treatment's
# sound_score object from Part 3a.
# One sig table per treatment x gene: DMR_boot_<TX>_<gene>_sig.

dmr_boots      <- list()
sig_list       <- list()
sig_treatments <- character(0)

for (i in seq_len(nrow(targets))) {
  trt  <- targets$Treatment[i]
  gene <- targets$Gene[i]
  key  <- paste0("DMR_boot_", trt, "_", gene)
  same_trt_genes <- targets$Gene[targets$Treatment == trt]

  dmr_boots[[key]] <- boot_score(sound_score_obj  = dmr_scores[[trt]],
                                 target_gene      = gene,
                                 scoring_col_name = scoring_col,
                                 target_start     = targets$target_start[i],
                                 target_end       = targets$target_end[i],
                                 exclude_genes    = if (exclude_other_targets)
                                                      setdiff(same_trt_genes, gene))
  sig_list[[paste0(key, "_sig")]] <- get_sig_dmrs(dmr_boots[[key]], min_score = min_sig_score)
  sig_treatments <- c(sig_treatments, trt)
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
