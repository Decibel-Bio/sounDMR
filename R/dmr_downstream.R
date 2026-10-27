#Copyright 2023 Sound Agriculture Company
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

# -----------------------------------------------------------------------------
# Downstream summaries and figures for significant DMRs.
#
# Replaces the ad-hoc combine / datasummary / DMR_Info / treatment_plots steps:
#   * nothing is selected by column position - everything by name
#   * treatment-specific column names (Z_TreatmentT3_small_mean, T3_mean_dif_mean,
#     ...) are mapped to shared names so any number of treatments rbind cleanly
#   * any number of target genes per treatment (summaries by Treatment and by
#     Treatment x Gene)
#   * every output column is atomic (no list columns)
#   * treatments / contexts with zero significant DMRs still appear in summaries
#
# Note: reshape2 is also imported by this package and masks data.table's
# dcast()/melt(), so those are always called as data.table::dcast / data.table::melt.
# -----------------------------------------------------------------------------

CX_LEVELS <- c("CG", "CHG", "CHH")
CX_COLORS <- c(CG = "#2a78d6", CHG = "#eb6834", CHH = "#1baf7a")

utils::globalVariables(c(
  ".", "..need", "..by", "..join_by", ".score", "CX", "Chr", "Count", "Cyt",
  "DMR_ID", "Gene", "Genome_High", "Genome_Low", "Group_Label", "Label", "Max",
  "Mid", "N", "Percent", "Percent_Control", "Percent_Target", "Percent_Dif_vs_Control",
  "Percent_Off_Target", "Percent_Dif_vs_Off_Target", "Off_Group",
  "Dif_Vs_Control_mean", "Dif_Vs_Off_Target_mean", "Off_Target_Fraction", "Off_Target_Adjusted_Score",
  "Percent_mean", "Plant", "Position", "RD_iqr", "RD_mean", "Role", "Score", "Source",
  "Start", "Stop", "Total_DMRs", "Total_Score", "Treatment", "cp_group", "i.DMR_ID",
  "score", "total_RD", "x.Group_Label", "x.Percent", "x.Plant", "x.Position",
  "x.total_RD", "xend", "y", "yend"
))


# -----------------------------------------------------------------------------
# 1. Selecting and combining significant DMR tables
# -----------------------------------------------------------------------------

#' Get Significant DMRs from a boot_score Object
#'
#' Pull the changepoint regions around the target gene whose precision-adjusted
#' score passes a threshold.
#'
#' @param boot_obj (list) - output of `boot_score()`
#' @param min_score (num) - keep regions with `score_col` strictly greater than this
#' @param score_col (str) - column to threshold on
#' @return (df) - the passing rows of `boot_obj$target_rs`, ungrouped
#' @examples
#' # DMR_boot_T2_sig <- get_sig_dmrs(DMR_boot_T2, min_score = 6)
#' @export
get_sig_dmrs <- function(boot_obj, min_score = 6, score_col = "adjusted_soundscore") {
  rs <- as.data.frame(boot_obj$target_rs)
  if (!score_col %in% names(rs)) stop("`", score_col, "` not found in boot_obj$target_rs")
  rs[!is.na(rs[[score_col]]) & rs[[score_col]] > min_score, , drop = FALSE]
}

#' @noRd
treatment_from_name <- function(list_name, name_pattern = "^DMR_boot_(.+)_sig$") {
  if (is.null(list_name) || !nzchar(list_name) || !grepl(name_pattern, list_name)) {
    return(NA_character_)
  }
  sub(name_pattern, "\\1", list_name)
}

#' @noRd
treatment_from_z <- function(df, target_variable = "Treatment", z_suffix = "_small") {
  pat <- paste0("^Z_", target_variable, "(.+)", z_suffix, "_mean$")
  hits <- grep(pat, names(df), value = TRUE)
  if (length(hits) == 1) sub(pat, "\\1", hits) else NA_character_
}

#' Infer the treatment a sig-DMR table belongs to
#'
#' Primary source is the object name (`DMR_boot_<TX>_sig` -> `<TX>`). The Z column
#' written by `sound_score()` (e.g. `Z_TreatmentT3_small_mean` -> `T3`) is used as
#' a check: if the name has extra text after the Z treatment
#' (`DMR_boot_T3_geneB_sig`) the Z treatment is used; if they disagree otherwise
#' the name is used with a warning; with no usable name the Z column is used.
#' @noRd
infer_treatment <- function(df, list_name = NULL, target_variable = "Treatment",
                            z_suffix = "_small", name_pattern = "^DMR_boot_(.+)_sig$") {
  from_name <- treatment_from_name(list_name, name_pattern)
  from_z    <- treatment_from_z(df, target_variable, z_suffix)

  if (is.na(from_name) && is.na(from_z)) {
    stop("Could not get a treatment for '", list_name, "': name doesn't match ",
         name_pattern, " and no Z_", target_variable, "..._mean column. ",
         "Pass `treatments =` explicitly.")
  }
  if (is.na(from_name)) return(from_z)
  if (is.na(from_z) || from_name == from_z) return(from_name)
  if (startsWith(from_name, paste0(from_z, "_"))) return(from_z)
  warning("'", list_name, "': name says treatment '", from_name,
          "' but its Z column is for '", from_z, "'. Using '", from_name, "'.")
  from_name
}

#' Map treatment-specific column names to shared ones
#'   MethRegion_Z_TreatmentT3_small_mean -> MethRegion_Z_mean
#'   T3_mean_dif_mean                    -> Dif_Vs_Control_mean
#'   T3_off_mean                         -> Dif_Vs_Off_Target_mean
#' @noRd
standardize_sig_cols <- function(dt, treatment, target_variable = "Treatment",
                                 z_suffix = "_small",
                                 per_change_suffix = "_mean_dif",
                                 off_suffix = "_off") {
  nm <- names(dt)
  nm <- gsub(paste0("Z_", target_variable, treatment, z_suffix), "Z", nm, fixed = TRUE)
  swap_prefix <- function(x, old, new) {
    hit <- startsWith(x, old)
    x[hit] <- paste0(new, substring(x[hit], nchar(old) + 1))
    x
  }
  nm <- swap_prefix(nm, paste0(treatment, per_change_suffix), "Dif_Vs_Control")
  nm <- swap_prefix(nm, paste0(treatment, off_suffix), "Dif_Vs_Off_Target")
  data.table::setnames(dt, nm)
  dt
}

#' Combine Significant DMR Tables
#'
#' Combine any number of significant-DMR tables (one per treatment, or one per
#' treatment x target gene) into one flat table with shared column names.
#'
#' The treatment for each table comes from its name (`DMR_boot_<TX>_sig` ->
#' `<TX>`), cross-checked against its `Z_<target_variable><TX>_small_mean` column.
#' Target gene(s) come from each table's `Gene` column, so a table can hold one or
#' several genes, and several tables can share a treatment.
#'
#' @param sig_list (named list of df) - e.g. `mget(ls(pattern = "^DMR_boot_.*_sig$"))`
#' @param treatments (chr) *optional* - treatment for each element of `sig_list`;
#' overrides inference
#' @param target_variable (str) - the fixed effect used in `find_DMR()`
#' @param z_suffix,per_change_suffix,off_suffix (str) - naming used for the Z,
#' percent-change and off-target columns
#' @param name_pattern (regex) - pattern with one capture group for the treatment
#' @param verbose (bool) - print the table -> treatment -> gene mapping
#' @return (data.table) - one row per DMR with `Treatment`, `Source` and `DMR_ID`
#' first. `attr(, "treatments")` holds every treatment seen (including ones with 0
#' DMRs), `attr(, "mapping")` the mapping that was used, and `attr(, "pairs")` every
#' treatment x gene, including empty tables whose gene is given by their name
#' (`DMR_boot_<TX>_<gene>_sig`).
#' @import data.table
#' @export
combine_sig_dmrs <- function(sig_list, treatments = NULL,
                             target_variable = "Treatment", z_suffix = "_small",
                             per_change_suffix = "_mean_dif", off_suffix = "_off",
                             name_pattern = "^DMR_boot_(.+)_sig$", verbose = TRUE) {
  if (is.data.frame(sig_list)) sig_list <- list(sig_list)
  if (!is.null(treatments) && length(treatments) != length(sig_list)) {
    stop("`treatments` must be the same length as `sig_list`.")
  }
  nms <- names(sig_list)
  if (is.null(nms)) nms <- rep("", length(sig_list))
  nms[is.na(nms)] <- ""
  src <- ifelse(nzchar(nms), nms, paste0("table_", seq_along(sig_list)))

  mapping <- vector("list", length(sig_list))
  pairs   <- vector("list", length(sig_list))
  pieces <- lapply(seq_along(sig_list), function(i) {
    df <- sig_list[[i]]
    trt <- if (!is.null(treatments)) treatments[i] else
      infer_treatment(df, nms[i], target_variable, z_suffix, name_pattern)
    genes <- if (nrow(df) && "Gene" %in% names(df)) sort(unique(as.character(df$Gene))) else character(0)
    # An empty table has no Gene column values; if it is named
    # DMR_boot_<TX>_<gene>_sig, take the gene from the name so the
    # treatment x gene pair still shows up (as zeros) in the by-gene summary.
    name_gene <- NA_character_
    if (!length(genes)) {
      stem <- paste0("DMR_boot_", trt, "_")
      if (startsWith(nms[i], stem) && endsWith(nms[i], "_sig")) {
        name_gene <- sub("_sig$", "", substring(nms[i], nchar(stem) + 1))
      }
    }
    mapping[[i]] <<- data.table(
      Source = src[i], Treatment = trt,
      Genes = if (length(genes)) paste(genes, collapse = "; ") else
        if (!is.na(name_gene)) paste(name_gene, "(no DMRs)") else "(no DMRs)",
      N_DMRs = nrow(df))
    pairs[[i]] <<- data.table(Treatment = trt,
                              Gene = if (length(genes)) genes else name_gene)
    if (nrow(df) == 0) return(NULL)
    dt <- as.data.table(as.data.frame(df))          # drops dplyr grouping
    # rename columns using the treatment they were actually written with
    col_trt <- treatment_from_z(df, target_variable, z_suffix)
    if (is.na(col_trt)) col_trt <- trt
    dt <- standardize_sig_cols(dt, col_trt, target_variable, z_suffix,
                               per_change_suffix, off_suffix)
    dt[, `:=`(Treatment = trt, Source = src[i])]
    dt
  })
  mapping <- rbindlist(mapping)
  if (verbose) {
    message("Sig table -> treatment -> gene mapping:")
    message(paste(utils::capture.output(print(mapping, row.names = FALSE)), collapse = "\n"))
  }

  out <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  if (nrow(out) == 0) stop("All sig-DMR tables are empty.")
  out[, DMR_ID := paste(Treatment, cp_group, sep = "|")]
  dup <- duplicated(out$DMR_ID)
  if (any(dup)) {
    message(sum(dup), " duplicate DMR rows (same treatment + cp_group) dropped.")
    out <- out[!dup]
  }
  setcolorder(out, c("Treatment", "Source", "DMR_ID"))
  setattr(out, "treatments", unique(mapping$Treatment))
  setattr(out, "mapping", mapping)
  pairs <- unique(rbindlist(pairs))
  setattr(out, "pairs", pairs[!is.na(Gene)])
  out[]
}


# -----------------------------------------------------------------------------
# 2. Per-DMR, per-treatment methylation and read depth (replaces DMR_Info)
# -----------------------------------------------------------------------------

#' Extract Cytosine-Level Data Inside Each DMR
#'
#' Pull every cytosine x sample observation in `LongPercent` that falls inside
#' each DMR (same gene, same context, between `Genome_Low` and `Genome_High`),
#' using a single non-equi join.
#'
#' @param dmrs (data.table) - output of `combine_sig_dmrs()`
#' @param long_percent (df) - `dmr_obj$LongPercent`
#' @param treatment_col (str) - column in `long_percent` holding treatment labels
#' @return (data.table) - `DMR_ID, Treatment, Plant, Position, Percent, total_RD`
#' @export
dmr_cytosine_data <- function(dmrs, long_percent, treatment_col = "Treatment") {
  need <- c("Gene", "CX", "Position", "Plant", "Percent", "total_RD", treatment_col)
  miss <- setdiff(need, names(long_percent))
  if (length(miss)) stop("long_percent is missing: ", paste(miss, collapse = ", "))

  lp <- as.data.table(long_percent)[Gene %in% unique(dmrs$Gene), ..need]
  setnames(lp, treatment_col, "Group_Label")
  regions <- unique(as.data.table(dmrs)[, .(DMR_ID, Gene, CX, Genome_Low, Genome_High)])

  lp[regions,
     on = .(Gene, CX, Position >= Genome_Low, Position <= Genome_High),
     nomatch = NULL, allow.cartesian = TRUE,
     .(DMR_ID = i.DMR_ID, Treatment = x.Group_Label, Plant = x.Plant,
       Position = x.Position, Percent = x.Percent, total_RD = x.total_RD)]
}

#' Per-DMR Statistics for Every Treatment Group
#'
#' @param cyto (data.table) - output of `dmr_cytosine_data()`
#' @return (data.table) - one row per DMR x treatment group: number of cytosines
#' and samples, mean and IQR of percent methylation and of read depth
#' @export
dmr_group_stats <- function(cyto) {
  cyto[, .(n_cytosines  = uniqueN(Position),
           n_samples    = uniqueN(Plant),
           Percent_mean = mean(Percent, na.rm = TRUE),
           Percent_iqr  = stats::IQR(Percent, na.rm = TRUE),
           RD_mean      = mean(total_RD, na.rm = TRUE),
           RD_iqr       = stats::IQR(total_RD, na.rm = TRUE)),
       by = .(DMR_ID, Treatment)]
}

#' Add Target, Control and Off-Target Methylation / Read Depth to the DMR Table
#'
#' Adds plain numeric columns measured inside each DMR:
#' \itemize{
#'   \item `Percent_Target`, `RD_Target`, `RD_Target_IQR` - the DMR's own treatment
#'   \item `Percent_Control`, `RD_Control`, `RD_Control_IQR` - the control group
#'   \item `Percent_Dif_vs_Control` - target minus control
#'   \item with `off_targets`: `Off_Targets` (which treatments), `Percent_Off_Target`,
#'   `RD_Off_Target` (mean over those treatment groups, each weighted equally) and
#'   `Percent_Dif_vs_Off_Target` - target minus off-target
#' }
#'
#' @param dmrs (data.table) - output of `combine_sig_dmrs()`
#' @param stats (data.table) - output of `dmr_group_stats()`
#' @param control (str) - control level of the treatment column
#' @param off_targets (named list) *optional* - for each treatment, the treatments
#' that targeted a DIFFERENT gene, used as an off-target baseline, e.g.
#' `list(T2 = c("T5", "T6", "T7"), T5 = c("T2", "T3", "T4"))`. Treatments not
#' named get no off-target columns (NA).
#' @return (data.table) - `dmrs` with the added columns
#' @export
add_target_stats <- function(dmrs, stats, control, off_targets = NULL) {
  if (!control %in% stats$Treatment) {
    warning("Control group '", control, "' not found in LongPercent; control columns will be NA.")
  }
  tgt <- stats[, .(DMR_ID, Treatment, Percent_Target = Percent_mean, RD_Target = RD_mean,
                   RD_Target_IQR = RD_iqr)]
  ctl <- stats[Treatment == control,
               .(DMR_ID, Percent_Control = Percent_mean, RD_Control = RD_mean,
                 RD_Control_IQR = RD_iqr)]
  out <- merge(as.data.table(dmrs), tgt, by = c("DMR_ID", "Treatment"), all.x = TRUE, sort = FALSE)
  out <- merge(out, ctl, by = "DMR_ID", all.x = TRUE, sort = FALSE)
  out[, Percent_Dif_vs_Control := Percent_Target - Percent_Control]

  ot <- off_target_table(off_targets)
  if (nrow(ot)) {
    missing_groups <- setdiff(ot$Off_Group, stats$Treatment)
    if (length(missing_groups)) {
      warning("Off-target treatment(s) not found in LongPercent: ",
              paste(missing_groups, collapse = ", "))
    }
    x <- merge(unique(out[, .(DMR_ID, Treatment)]), ot, by = "Treatment", allow.cartesian = TRUE)
    x <- merge(x, stats[, .(DMR_ID, Off_Group = Treatment, Percent_mean, RD_mean)],
               by = c("DMR_ID", "Off_Group"))
    offs <- x[, .(Off_Targets        = paste(sort(unique(Off_Group)), collapse = "; "),
                  Percent_Off_Target = mean(Percent_mean, na.rm = TRUE),
                  RD_Off_Target      = mean(RD_mean, na.rm = TRUE)),
              by = DMR_ID]
    out <- merge(out, offs, by = "DMR_ID", all.x = TRUE, sort = FALSE)
    out[, Percent_Dif_vs_Off_Target := Percent_Target - Percent_Off_Target]
  }

  if (anyNA(out$RD_Target)) {
    warning(sum(is.na(out$RD_Target)), " DMR(s) had no matching rows in LongPercent ",
            "- check that treatment labels match.")
  }
  setcolorder(out, names(dmrs))
  out[]
}

#' Off-Target Fraction and Off-Target Adjusted Score
#'
#' For DMRs whose treatment has off-target treatments (so the region summary has
#' `Dif_Vs_Off_Target_mean`), adds:
#' \itemize{
#'   \item `Off_Target_Fraction` = (Dif_Vs_Control_mean - Dif_Vs_Off_Target_mean) /
#'   Dif_Vs_Control_mean, which equals (off-target - control) / (treatment - control):
#'   the share of the treatment's change that the off-target treatments also show.
#'   0 = fully specific to the targeting treatment; 1 = off-targets changed just as much;
#'   below 0 = off-targets moved the other way.
#'   \item `Off_Target_Adjusted_Score` = score x (1 + (1 - Off_Target_Fraction)):
#'   x2 for a fully specific DMR, x1 for a non-specific one.
#' }
#' DMRs without off-targets get NA in both.
#'
#' @param dmrs (data.table) - combined DMR table (from `combine_sig_dmrs()`)
#' @param score_col (str) - score to adjust, usually `adjusted_soundscore`
#' @return (data.table) - `dmrs` with the two columns added (if
#' `Dif_Vs_Off_Target_mean` is present; otherwise unchanged)
#' @export
add_off_target_scores <- function(dmrs, score_col = "adjusted_soundscore") {
  out <- as.data.table(dmrs)
  if (!all(c("Dif_Vs_Control_mean", "Dif_Vs_Off_Target_mean") %in% names(out))) return(out[])
  out[, Off_Target_Fraction := ifelse(is.na(Dif_Vs_Control_mean) | Dif_Vs_Control_mean == 0,
                                      NA_real_,
                                      (Dif_Vs_Control_mean - Dif_Vs_Off_Target_mean) / Dif_Vs_Control_mean)]
  out[, Off_Target_Adjusted_Score := get(score_col) * (1 + (1 - Off_Target_Fraction))]
  out[]
}

#' Long table of treatment -> off-target treatment pairs (self-references dropped)
#' @noRd
off_target_table <- function(off_targets) {
  if (is.null(off_targets) || !length(off_targets)) {
    return(data.table(Treatment = character(0), Off_Group = character(0)))
  }
  if (is.null(names(off_targets)) || any(!nzchar(names(off_targets)))) {
    stop("`off_targets` must be a named list: list(<treatment> = c(<off-target treatments>), ...)")
  }
  ot <- rbindlist(lapply(names(off_targets), function(t) {
    g <- setdiff(as.character(off_targets[[t]]), t)
    if (length(g)) data.table(Treatment = t, Off_Group = g) else NULL
  }))
  if (!nrow(ot)) ot <- data.table(Treatment = character(0), Off_Group = character(0))
  unique(ot)
}


# -----------------------------------------------------------------------------
# 3. Summary tables
# -----------------------------------------------------------------------------

#' Summarise DMR Scores by Treatment and Context
#'
#' Sum, count and max of the DMR score per cytosine context, for each treatment
#' (or treatment x gene). Every context always gets a column and every treatment in
#' `treatments` gets a row, even with no significant DMRs.
#'
#' @param dmrs (data.table) - combined DMR table
#' @param by (chr) - grouping columns, e.g. `"Treatment"` or `c("Treatment", "Gene")`
#' @param score_col (str) - score to sum / max
#' @param treatments (chr) - every treatment that should appear
#' @param cx_levels (chr) - cytosine contexts
#' @param design (df) *optional* - `Treatment` (and optionally `Gene`) plus
#' annotation columns (e.g. `Description`, `Oligo_Duplexes`) to merge in
#' @param required_keys (df) *optional* - combinations of the `by` columns that must
#' appear even with no significant DMRs, e.g. `attr(dmrs, "pairs")` for
#' treatment x gene summaries
#' @return (data.table) - summary table
#' @export
summarise_dmrs <- function(dmrs, by = "Treatment", score_col = "adjusted_soundscore",
                           treatments = attr(dmrs, "treatments"),
                           cx_levels = CX_LEVELS, design = NULL, required_keys = NULL) {
  dt <- as.data.table(dmrs)
  dt[, .score := get(score_col)]

  long <- dt[, .(Score = sum(.score, na.rm = TRUE), N = .N,
                 Max = max(.score, na.rm = TRUE), Cyt = sum(Count, na.rm = TRUE)),
             by = c(by, "CX")]

  # full grid so empty contexts / treatments are explicit zeros, not missing columns
  keys <- unique(dt[, ..by])
  if (!is.null(treatments) && "Treatment" %in% by) {
    extra <- setdiff(treatments, keys$Treatment)
    if (length(extra)) keys <- rbind(keys, data.table(Treatment = extra), fill = TRUE)
  }
  if (!is.null(required_keys) && all(by %in% names(required_keys))) {
    keys <- unique(rbind(keys, unique(as.data.table(required_keys)[, by, with = FALSE])))
  }
  grid <- keys[, .(CX = cx_levels), by = by]
  long <- merge(grid, long, by = c(by, "CX"), all.x = TRUE)
  long[is.na(Score), `:=`(Score = 0, N = 0L, Cyt = 0L)]
  long[, CX := factor(CX, levels = cx_levels)]

  wide <- data.table::dcast(long, stats::as.formula(paste(paste(by, collapse = " + "), "~ CX")),
                            value.var = c("Score", "N", "Max"), sep = "_")
  totals <- long[, .(Total_Score = sum(Score), Total_DMRs = sum(N),
                     Cytosines_Covered = sum(Cyt)), by = by]
  out <- merge(wide, totals, by = by)

  # per-treatment summaries: list the target gene(s) each treatment's DMRs came from
  if (!"Gene" %in% by && "Gene" %in% names(dt)) {
    genes <- dt[, .(Genes = paste(sort(unique(Gene)), collapse = "; ")), by = by]
    out <- merge(out, genes, by = by, all.x = TRUE)
  }

  if (!is.null(design)) {
    design <- as.data.table(design)
    join_by <- intersect(by, names(design))
    if (!"Treatment" %in% join_by) stop("`design` needs a Treatment column.")
    unmatched <- setdiff(design$Treatment, out$Treatment)
    if (length(unmatched)) {
      warning("design Treatment(s) not found in the DMR tables: ",
              paste(unmatched, collapse = ", "),
              " - check spelling against LongPercent$Treatment.")
    }
    # design may be finer than this summary (e.g. per gene while summarising per
    # treatment): collapse to one row per key so rows are never duplicated
    extra_keys <- setdiff(intersect(names(design), c("Treatment", "Gene")), join_by)
    if (length(extra_keys) || anyDuplicated(design[, ..join_by])) {
      ann <- setdiff(names(design), c(join_by, extra_keys))
      design <- design[, lapply(.SD, function(x) {
        x <- unique(stats::na.omit(as.character(x)))
        if (length(x)) paste(x, collapse = "; ") else NA_character_
      }), by = join_by, .SDcols = ann]
    }
    out <- merge(out, design, by = join_by, all.x = TRUE)
  } else {
    out[, `:=`(Description = NA_character_, Oligo_Duplexes = NA_character_)]
  }

  col_order <- c(by, "Genes",
                 paste0("Score_", cx_levels), paste0("N_", cx_levels),
                 paste0("Max_", cx_levels),
                 "Total_Score", "Total_DMRs", "Cytosines_Covered")
  setcolorder(out, intersect(col_order, names(out)))
  num <- names(out)[vapply(out, is.numeric, logical(1))]
  out[, (num) := lapply(.SD, round, 2), .SDcols = num]
  setorderv(out, by)
  out[]
}

#' Human-Readable Summary Column Names
#'
#' `Score_CG` -> `CG Score`, `N_CG` -> `# CG DMRs`, `Max_CG` -> `Max CG Score`,
#' underscores -> spaces.
#'
#' @param summ (data.table) - output of `summarise_dmrs()`
#' @return (data.table) - renamed copy
#' @export
prettify_summary <- function(summ) {
  s <- copy(as.data.table(summ))
  nm <- names(s)
  nm <- sub("^Score_(.+)$", "\\1 Score", nm)
  nm <- sub("^N_(.+)$", "# \\1 DMRs", nm)
  nm <- sub("^Max_(.+)$", "Max \\1 Score", nm)
  nm <- gsub("_", " ", nm)
  setnames(s, nm)
  s[]
}


# -----------------------------------------------------------------------------
# 4. BED export
# -----------------------------------------------------------------------------

#' Write Significant DMRs as a BED File
#'
#' Writes a valid BED4 (tab-separated, no header, 0-based half-open). Positions in
#' sounDMR come from methyl-bed column 2 (0-based start), so by default
#' start = `Genome_Low` and end = `Genome_High + 1`.
#'
#' @param dmrs (data.table) - combined DMR table
#' @param file (str) - output path
#' @param score_col (str) - score included in the feature name
#' @param positions_0based (bool) - set FALSE if `Position` is 1-based
#' @return (data.table) - the BED table, invisibly
#' @export
write_dmr_bed <- function(dmrs, file, score_col = "adjusted_soundscore",
                          positions_0based = TRUE) {
  dt <- as.data.table(dmrs)
  multi_gene <- dt[, uniqueN(Gene) > 1]
  label <- if (multi_gene) paste(dt$Treatment, dt$Gene, dt$CX, sep = "_") else
    paste(dt$Treatment, dt$CX, sep = "_")
  bed <- data.table(
    chrom = dt$Chr,
    start = if (positions_0based) dt$Genome_Low else dt$Genome_Low - 1L,
    end   = if (positions_0based) dt$Genome_High + 1L else dt$Genome_High,
    name  = paste(label, round(dt[[score_col]], 2), sep = "_")
  )
  setorderv(bed, c("chrom", "start"))
  fwrite(bed, file, sep = "\t", col.names = FALSE, quote = FALSE)
  invisible(bed)
}


# -----------------------------------------------------------------------------
# 5. Figures
# -----------------------------------------------------------------------------

#' @noRd
theme_dmr <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(colour = "grey35", size = 9),
          axis.text = element_text(colour = "grey30"),
          strip.text = element_text(face = "bold"))
}

#' Plot Methylation in One DMR Across Treatment Groups
#'
#' Boxes are per-cytosine group means; dots are each plant's mean over the DMR.
#' The DMR's own treatment is blue, the control grey, and its off-target
#' treatments (the `Off_Targets` column added by `add_target_stats()`) orange.
#' Replaces `cpt_compare()` / `treatment_plots()`.
#'
#' @param dmr_id (str) - a `DMR_ID` from the combined table
#' @param cyto (data.table) - output of `dmr_cytosine_data()`
#' @param dmrs (data.table) - combined DMR table
#' @param control (str) - control level
#' @param treatment_levels (chr) *optional* - x-axis order
#' @return ggplot object
#' @import ggplot2
#' @export
plot_dmr_methylation <- function(dmr_id, cyto, dmrs, control, treatment_levels = NULL) {
  d <- cyto[DMR_ID == dmr_id]
  info <- as.data.table(dmrs)[DMR_ID == dmr_id][1]
  if (is.null(treatment_levels)) {
    treatment_levels <- c(control, sort(setdiff(unique(d$Treatment), control)))
  }
  d[, Treatment := factor(Treatment, levels = treatment_levels)]
  per_cyt   <- d[, .(Percent = mean(Percent, na.rm = TRUE)), by = .(Treatment, Position)]
  per_plant <- d[, .(Percent = mean(Percent, na.rm = TRUE)), by = .(Treatment, Plant)]
  off_trts <- if ("Off_Targets" %in% names(info) && !is.na(info$Off_Targets)) {
    strsplit(info$Off_Targets, "; ", fixed = TRUE)[[1]]
  } else character(0)
  role_levels <- c("Target treatment", "Control", "Off-target treatments", "Other treatments")
  role <- function(t) factor(ifelse(t == info$Treatment, role_levels[1],
                             ifelse(t == control, role_levels[2],
                             ifelse(t %in% off_trts, role_levels[3], role_levels[4]))),
                             levels = role_levels)
  per_cyt[, Role := role(Treatment)]
  per_plant[, Role := role(Treatment)]
  fills <- stats::setNames(c("#2a78d6", "#8c8c8c", "#f2a07b", "#e3e3e3"), role_levels)

  ggplot(per_cyt, aes(Treatment, Percent)) +
    geom_boxplot(aes(fill = Role), outlier.shape = NA, width = 0.6,
                 colour = "grey30", linewidth = 0.3, na.rm = TRUE) +
    geom_point(data = per_plant, size = 2.2, shape = 21, fill = "white",
               colour = "grey15", stroke = 0.6, na.rm = TRUE) +
    scale_fill_manual(values = fills, drop = TRUE, name = NULL) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
    labs(title = paste0(info$Treatment, "  ", info$cp_group),
         subtitle = sprintf("%s  %s:%s-%s  |  %d cytosines  |  adj. score %.2f  |  box = per-cytosine group mean, dots = per-plant mean",
                            info$CX, info$Chr, format(info$Genome_Low, big.mark = ","),
                            format(info$Genome_High, big.mark = ","), info$Count,
                            info$adjusted_soundscore),
         x = NULL, y = "% methylation") +
    theme_dmr() +
    theme(legend.position = "bottom")
}

#' Map of Significant DMRs Relative to the Target
#'
#' One panel per gene x context; line = DMR extent, point size = score; dashed
#' line at position 0.
#'
#' @param dmrs (data.table) - combined DMR table
#' @param score_col (str) - score used for point size
#' @param treatment_levels (chr) *optional* - y-axis order
#' @return ggplot object
#' @export
plot_dmr_map <- function(dmrs, score_col = "adjusted_soundscore",
                         treatment_levels = NULL) {
  d <- copy(as.data.table(dmrs))
  d[, score := get(score_col)]
  d[, Mid := (Start + Stop) / 2]
  if (is.null(treatment_levels)) treatment_levels <- sort(unique(d$Treatment))
  d[, Treatment := factor(Treatment, levels = rev(treatment_levels))]
  d[, CX := factor(CX, levels = CX_LEVELS)]

  ggplot(d, aes(y = Treatment, colour = CX)) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_segment(aes(x = Start, xend = Stop, yend = Treatment), linewidth = 2,
                 lineend = "round", alpha = 0.6) +
    geom_point(aes(x = Mid, size = score), alpha = 0.9) +
    scale_colour_manual(values = CX_COLORS, drop = FALSE, guide = "none") +
    scale_size_area(max_size = 6, name = "Adj. score") +
    facet_grid(CX ~ Gene, scales = "free_x", drop = FALSE) +
    labs(x = "Position relative to target (bp)", y = NULL,
         title = "Significant DMRs by treatment") +
    theme_dmr() +
    theme(panel.grid.major.y = element_line(colour = "grey92"),
          panel.border = element_rect(fill = NA, colour = "grey85"))
}

#' Bar Chart of Summed DMR Score per Treatment
#'
#' @param summ (data.table) - output of `summarise_dmrs()`
#' @param by (chr) - the grouping used to build `summ`
#' @return ggplot object
#' @export
plot_score_summary <- function(summ, by = "Treatment") {
  s <- as.data.table(summ)
  long <- data.table::melt(s, id.vars = by, measure.vars = paste0("Score_", CX_LEVELS),
                           variable.name = "CX", value.name = "Score")
  long[, CX := factor(sub("^Score_", "", CX), levels = rev(CX_LEVELS))]
  tot <- s[, c(by, "Total_Score", "Total_DMRs"), with = FALSE]
  if (length(by) > 1) {
    long[, Label := do.call(paste, c(.SD, sep = "\n")), .SDcols = by]
    tot[, Label := do.call(paste, c(.SD, sep = "\n")), .SDcols = by]
  } else {
    long[, Label := get(by)]
    tot[, Label := get(by)]
  }
  ggplot(long, aes(Label, Score, fill = CX)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.5) +
    geom_text(data = tot, aes(Label, Total_Score,
                              label = sprintf("%.1f  (%d DMRs)", Total_Score, Total_DMRs)),
              inherit.aes = FALSE, vjust = -0.5, size = 3, colour = "grey25") +
    scale_fill_manual(values = CX_COLORS, breaks = CX_LEVELS, name = "Context") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Summed adjusted score", title = "DMR score by treatment") +
    theme_dmr()
}


# -----------------------------------------------------------------------------
# 6. One-call wrapper
# -----------------------------------------------------------------------------

#' Run the Downstream DMR Summary Pipeline
#'
#' Combines significant-DMR tables, adds per-treatment methylation and read depth
#' for every DMR, builds summary tables, and writes CSVs, a BED file and figures.
#' Supersedes the manual combine / `datasummary()` / `DMR_Info()` /
#' `treatment_plots()` steps.
#'
#' Files written to `out_dir`, each starting with `prefix`:
#' \itemize{
#'   \item `_sig_table_mapping.csv` - which sig table -> treatment -> gene(s)
#'   \item `_dmrs.csv` - one row per DMR, with target / control (/ off-target) read
#'   depth and \% methylation
#'   \item `_dmr_group_stats.csv` - one row per DMR x treatment group
#'   \item `_summary.csv` - per-treatment summary
#'   \item `_summary_by_gene.csv` - per-treatment x gene summary
#'   \item `_dmrs.bed` - BED4 of significant DMRs
#'   \item `_score_summary.pdf`, `_dmr_map.pdf`, `_dmr_methylation.pdf` (one page per DMR)
#' }
#'
#' @param sig_list (named list of df) - significant DMR tables, see `combine_sig_dmrs()`
#' @param dmr_obj (list) - output of `create_dmr_obj()`
#' @param control (str) - control level of `target_variable` in `dmr_obj$LongPercent`
#' @param out_dir (str) - output directory (created if needed)
#' @param prefix (str) - file name prefix, e.g. the experiment ID
#' @param treatments (chr) *optional* - treatment for each element of `sig_list`
#' @param design (df) *optional* - annotation merged into the summaries, see
#' `summarise_dmrs()`
#' @param off_targets (named list) *optional* - for each treatment, the treatments
#' that targeted a different gene (off-target baseline), e.g.
#' `list(T2 = c("T5", "T6", "T7"))`. Adds off-target % methylation and
#' difference columns to `_dmrs.csv`, plus `Off_Target_Fraction` and
#' `Off_Target_Adjusted_Score` (see `add_off_target_scores()`), and colours the
#' off-target treatments in the methylation plots.
#' See `add_target_stats()`.
#' @param target_variable (str) - the fixed effect used in `find_DMR()`
#' @param score_col (str) - DMR score column
#' @param positions_0based (bool) - see `write_dmr_bed()`
#' @param make_plots (bool) - write the PDF figures
#' @param per_dmr_pdfs (bool) - also write one PDF per DMR
#' @return (list), invisibly: `mapping`, `dmrs`, `group_stats`, `cytosines`,
#' `summary`, `summary_by_gene`
#' @examples
#' # sig_list <- mget(ls(pattern = "^DMR_boot_.*_sig$"))
#' # res <- run_dmr_downstream(sig_list, dmr_obj, control = "Water",
#' #                           out_dir = "db210_downstream", prefix = "db210")
#' @export
run_dmr_downstream <- function(sig_list, dmr_obj, control, out_dir = ".", prefix = "dmr",
                               treatments = NULL, design = NULL, off_targets = NULL,
                               target_variable = "Treatment",
                               score_col = "adjusted_soundscore",
                               positions_0based = TRUE, make_plots = TRUE,
                               per_dmr_pdfs = FALSE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  f <- function(x) file.path(out_dir, paste0(prefix, x))

  message("Combining sig-DMR tables")
  dmrs <- combine_sig_dmrs(sig_list, treatments, target_variable)
  all_trt <- attr(dmrs, "treatments")
  mapping <- attr(dmrs, "mapping")
  pairs   <- attr(dmrs, "pairs")

  message("Extracting per-cytosine data for ", nrow(dmrs), " DMRs")
  cyto  <- dmr_cytosine_data(dmrs, dmr_obj$LongPercent, treatment_col = target_variable)
  stats <- dmr_group_stats(cyto)
  dmrs  <- add_target_stats(dmrs, stats, control, off_targets)
  dmrs  <- add_off_target_scores(dmrs, score_col)   # no-op when no off-targets
  setattr(dmrs, "treatments", all_trt)
  setattr(dmrs, "mapping", mapping)

  message("Building summaries")
  summ_trt  <- summarise_dmrs(dmrs, by = "Treatment", score_col = score_col,
                              treatments = all_trt, design = design)
  summ_gene <- summarise_dmrs(dmrs, by = c("Treatment", "Gene"), score_col = score_col,
                              treatments = NULL, design = design, required_keys = pairs)

  fwrite(mapping, f("_sig_table_mapping.csv"))
  fwrite(dmrs, f("_dmrs.csv"))
  fwrite(stats, f("_dmr_group_stats.csv"))
  fwrite(prettify_summary(summ_trt), f("_summary.csv"))
  fwrite(prettify_summary(summ_gene), f("_summary_by_gene.csv"))
  write_dmr_bed(dmrs, f("_dmrs.bed"), score_col, positions_0based)

  if (make_plots) {
    message("Drawing figures")
    n_genes <- uniqueN(dmrs$Gene)
    ggsave(f("_score_summary.pdf"), plot_score_summary(summ_trt),
           width = max(5, 0.9 * nrow(summ_trt) + 2), height = 4.5)
    ggsave(f("_dmr_map.pdf"), plot_dmr_map(dmrs, score_col),
           width = max(7, 4 * n_genes + 2), height = 2 + 0.35 * length(all_trt) * 3)

    trt_levels <- c(control, sort(setdiff(unique(cyto$Treatment), control)))
    ord <- dmrs[order(Treatment, Gene, CX, Start), DMR_ID]
    grDevices::pdf(f("_dmr_methylation.pdf"), width = 8, height = 4.5)
    for (id in ord) print(plot_dmr_methylation(id, cyto, dmrs, control, trt_levels))
    grDevices::dev.off()
    if (per_dmr_pdfs) {
      for (id in ord) {
        ggsave(f(paste0("_", gsub("[|:]", "_", id), ".pdf")),
               plot_dmr_methylation(id, cyto, dmrs, control, trt_levels),
               width = 8, height = 4.5)
      }
    }
  }

  message("Done: outputs in ", normalizePath(out_dir))
  invisible(list(mapping = mapping, dmrs = dmrs, group_stats = stats, cytosines = cyto,
                 summary = summ_trt, summary_by_gene = summ_gene))
}
