# sounDMR example workflows

Pick the script that matches how your treatments map to target genes.

| Case | Design | Example | Script |
|---|---|---|---|
| **1 (default)** | One target gene, shared by all treatments | db210: T2–T5 all target EAB_chr6B | `DMR_workflow.R` |
| **2** | One treatment targets several genes | db187: T2_PAL2 targets LsPAL1 + LsPAL2 | `DMR_workflow_multi_gene.R` |
| **3** | Different treatments target different genes (a treatment may still have several) | T2 → LsPAL1, T3 → LsPAL2 | `DMR_workflow_treatment_specific_genes.R` |

**Step-through version (any case): `DMR_workflow_stepwise.R`.** Instead of looping
over treatments, Part 3 is a block you run for one treatment at a time. You set that
treatment's changepoint penalties, target gene(s) and window, score direction and
significance cutoff, check the plots and tables, rerun if needed, then move on to the
next treatment. This is closest to the original hand-run scripts. It handles all three
cases through the `targets_trt` table in step 3.1: one row for cases 1 and 3, several
rows for case 2. The settings used for each treatment are saved as
`<id>_settings_log.csv`. It comes set up for db210 (T2–T5 at EAB_chr6B, T6–T9 as
off-targets), like `DMR_workflow.R`.

`DMR_analysis_workflow.R` is the original generic example (Group C vs T, whole-genome
option, individual DMR model).

## What changes between cases

| Step | Case 1: one gene | Case 2: one treatment, several genes | Case 3: gene(s) per treatment |
|---|---|---|---|
| **Settings** | `target_gene`, `target_start`, `target_end` | `treatments_to_score` + a `target_genes` table (Gene, window) | a `targets` table (Treatment, Gene, window), one row per pair |
| **Gene coordinates (Geneco)** | Target gene plus background genes | All target genes plus background genes | All target genes plus background genes |
| **Part 1–2**: load data, `create_dmr_obj()`, `create_methyl_summary()`, `find_DMR()` | Same | Same | Same |
| **Part 3a**: `changepoint_analysis()` + `sound_score()` | Once per treatment | Once per treatment, **reused for all of its genes** | Once per treatment |
| **Part 3b**: `boot_score()` | Once per treatment, all at the same gene | Once per treatment × gene; `exclude_genes` = the treatment's **other** target genes | Once per `targets` row; `exclude_genes` = that treatment's other targets (if any) |
| **Sig table names** | `DMR_boot_<TX>_sig` | `DMR_boot_<TX>_<gene>_sig` | `DMR_boot_<TX>_<gene>_sig` |
| **Part 4.0**: final significance cutoff, one for all or per treatment (`new_min_sig_score`, default 8; `NULL` keeps `min_sig_score`) | Same | Same | Same |
| **Part 4**: `run_dmr_downstream()` | Same call | Same call | Same call |
| **Main summary table** | `<id>_summary.csv` (one row per treatment) | `<id>_summary_by_gene.csv` | `<id>_summary_by_gene.csv` |
| **Stepwise script, step 3.1** | Same one-row `targets_trt` in every block | Several rows in that treatment's block | A different one-row `targets_trt` per block |
| **Off-target treatments** (optional, see below) | `off_target_trts`: one set, used for every scored treatment | `off_target_trts`: one set, used for every scored treatment | `off_targets`: built per treatment from `targets` (editable) |
| **Off-target treatments, stepwise script** | `off_target_trts` in step 3.1 of each block | `off_target_trts` in step 3.1 of each block | `off_target_trts` in step 3.1 of each block |

## Off-target treatments

When different treatments target different genes, the treatments that hit *another*
gene are a useful second baseline. In addition to treated vs. control, you can ask
whether the change at this gene is specific to the treatments that targeted it.

Example: T2, T3 and T4 targeted gene A, and T5, T6 and T7 targeted gene B. When scoring
T2 at gene A, its off-target treatments are T5, T6 and T7, **not** T3 or T4 (they hit
gene A too).

| Script | Where you set it | Example |
|---|---|---|
| `DMR_workflow.R`, `DMR_workflow_multi_gene.R` | `off_target_trts` in SETTINGS (same set for every scored treatment) | `treatments_to_score <- c("T2","T3","T4")`, `off_target_trts <- c("T5","T6","T7")` |
| `DMR_workflow_treatment_specific_genes.R` | `off_targets` in SETTINGS, built automatically from `targets` (each treatment → the treatments sharing none of its genes); edit any entry by hand | `T2 → T5, T6, T7`; `T5 → T2, T3, T4` |
| `DMR_workflow_stepwise.R` | `off_target_trts` in step 3.1, chosen per treatment and recorded in `_settings_log.csv` | db210, running T2–T5: `off_target_trts <- c("T6","T7","T8","T9")` |

What you get, alongside every "vs. control" value:

| Level | vs. control | vs. off-target treatments |
|---|---|---|
| Per cytosine (`methyl_summary`) | `<TX>_mean_dif` = `<TX>_mean` − control mean | `<TX>_off` = `<TX>_mean` − mean of the off-target treatments' means |
| Per region (`sound_score()` region summary) | `<TX>_mean_dif_mean` / `_max` | `<TX>_off_mean` / `_max` |
| Per significant DMR (`<id>_dmrs.csv`) | `Dif_Vs_Control_mean`; observed `Percent_Control`, `Percent_Dif_vs_Control` | `Dif_Vs_Off_Target_mean`; observed `Off_Targets`, `Percent_Off_Target`, `RD_Off_Target`, `Percent_Dif_vs_Off_Target`; `Off_Target_Fraction` = (Dif_Vs_Control − Dif_Vs_Off_Target) / Dif_Vs_Control and `Off_Target_Adjusted_Score` = score × (1 + (1 − Off_Target_Fraction)) |
| Figures (`<id>_dmr_methylation.pdf`) | control box in grey | off-target boxes in orange |

The off-target treatments are averaged with equal weight per treatment. The DMR scores
and the significance cutoff are still based on the treatment vs. control model; the
off-target columns are there to compare against, not to select DMRs.

## Why `exclude_genes` matters in case 2

`boot_score()` compares the target region's score against random positions in
the *other* genes. If the same treatment also targeted some of those genes,
their real DMRs end up in the background, which makes the bootstrap
conservative. `exclude_genes` (new in 3.1.1) drops them from the background.
In case 3, genes targeted only by *other* treatments stay in the background,
since this treatment should have no effect there.

## Outputs (all cases)

`run_dmr_downstream()` writes these to `out_dir`, each prefixed with the experiment ID:

- `_sig_table_mapping.csv`: which sig table → treatment → gene(s)
- `_dmrs.csv`: one row per significant DMR, with target/control % methylation and read depth
- `_dmr_group_stats.csv`: one row per DMR × treatment group
- `_summary.csv`, `_summary_by_gene.csv`: scores, DMR counts and max score by context
- `_dmrs.bed`: BED4 of the significant DMRs
- `_score_summary.pdf`, `_dmr_map.pdf`, `_dmr_methylation.pdf`: figures

At the end, each example script also writes:

- `<id>_ms.csv`: the combined `methyl_summary`, with every treatment's columns
- `<id>.RData`: the whole session (`save.image()`); reload with `load()`
