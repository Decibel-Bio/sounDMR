# sounDMR 3.1.2
	- Off-target treatments: run_dmr_downstream() / add_target_stats() take `off_targets`,
	  a named list of the treatments that targeted a DIFFERENT gene for each treatment
	  (e.g. list(T2 = c("T5","T6","T7"))). Adds Off_Targets, Percent_Off_Target,
	  RD_Off_Target and Percent_Dif_vs_Off_Target to _dmrs.csv; plot_dmr_methylation()
	  shows off-target treatments in orange with a legend
	- sound_score(): new dmr_score4 (+ dmr_score4_Percentile) = dmr_score3 with region size
	  n^(1/2.5) instead of n^(1/3), giving larger regions a slight boost; exponent set by
	  `score4_count_exponent`. Signal gate: dmr_score4 is 0 for regions with |mean Z| <
	  `score4_min_mean_z` (0.25) or |mean change| < `score4_min_mean_dif` (2 % points), so
	  long flat regions can't score on size alone (flag in dmr_score4_signal). Existing
	  scores are unchanged
	- add_off_target_scores() (run automatically by run_dmr_downstream() when off-targets
	  are set): Off_Target_Fraction = (Dif_Vs_Control - Dif_Vs_Off_Target) / Dif_Vs_Control
	  and Off_Target_Adjusted_Score = score x (1 + (1 - Off_Target_Fraction)) in _dmrs.csv
	- Renamed for clarity: Percent_Dif_Observed -> Percent_Dif_vs_Control,
	  Off_Target_Dif_* -> Dif_Vs_Off_Target_*
	- Example workflows: `use_off_target` (all other treatments) replaced by explicit
	  off-target treatment sets (off_target_trts / off_targets); see README.md
	- DMR_workflow_stepwise.R: all treatments share one methyl_summary; optional step 4.0
	  re-applies a significance cutoff (one for all, or per treatment) to the kept
	  treatments without rerunning Part 3; step 3.1 check no longer fails when
	  off_target_trts is NULL
	- All example workflows (not just stepwise) have an optional step 4.0 to re-apply a
	  significance cutoff, one for all or per treatment, before run_dmr_downstream()
	- DMR_workflow.R (default) and DMR_workflow_stepwise.R are set up for db210: T2-T5 at
	  EAB_chr6B, T6-T9 as the off-target control; reproduces the saved db210 sig tables
	- Example workflows end by writing the combined methyl_summary to <id>_ms.csv and
	  saving the whole session to <id>.RData (replaces the _results.rds file)

# sounDMR 3.1.1
	- boot_score(): new optional `exclude_genes` argument to leave other genes targeted
	  by the same treatment out of the bootstrap background (default unchanged)
	- Example workflows split by experimental design (inst/examples, see README.md):
		- DMR_workflow.R: one target gene shared by all treatments (default)
		- DMR_workflow_multi_gene.R: one treatment targets several genes
		- DMR_workflow_treatment_specific_genes.R: treatments target different genes
		- DMR_workflow_stepwise.R: any design, run one treatment at a time with
		  treatment-specific penalties, windows, direction and cutoffs; writes a settings log
	- Replaces inst/examples/DMR_multi_target_workflow.R from 3.1.0
	- boot_score(): when several regions tie for the target's max score, take the first
	  (as the bootstrap rows already did) instead of warning
	- run_dmr_downstream(): a treatment x gene with no significant DMRs now appears as a
	  zero row in _summary_by_gene.csv (gene taken from the table name
	  DMR_boot_<TX>_<gene>_sig); summarise_dmrs() gains `required_keys`

# sounDMR 3.1.0
	- New downstream pipeline for significant DMRs (R/dmr_downstream.R):
		- run_dmr_downstream(): one call that combines sig-DMR tables and writes summary
		  CSVs, a BED file and figures
		- get_sig_dmrs(): threshold a boot_score() object's target regions
		- combine_sig_dmrs(): combines any number of sig tables by column name; the
		  treatment comes from the object name (DMR_boot_<TX>_sig) and the gene(s) from
		  the Gene column, so several target genes per treatment are supported
		- dmr_cytosine_data(), dmr_group_stats(), add_target_stats(): per-DMR, per-group
		  methylation and read depth as plain numeric columns (fixes the empty list-column
		  RD_Target produced by DMR_Info())
		- summarise_dmrs(), prettify_summary(): summary tables by treatment and by
		  treatment x gene; empty contexts / treatments appear as zeros
		- write_dmr_bed(): valid BED4 output (tab-separated, 0-based half-open)
		- plot_dmr_methylation(), plot_dmr_map(), plot_score_summary(): figures
	- DMR_Info() and treatment_plots() are kept for compatibility but superseded
	- NAMESPACE regenerated with roxygen2; man pages added for Treatment_Compare(),
	  cpt_compare(), treatment_plots() and DMR_Info()

# sounDMR 2.3.2
	- Added new feature to calculate number of cytosines included in a changepoint region
	- Fixed a bug in creating line plots for changepoints

# sounDMR 2.3.1
        - In create_dmr_obj():
              - Changed calculation of Percent column in LongPercent to weighted mean
              - Fixed bug in LongMeth total_RD calculation to use sum of read depths instead of mean

# sounDMR 2.3.0
	- Removed split_by_chromosome function from workflow due to bugs 
	- Improved notes for running whole genome analysis
	- Improved memory issues with group_DMR()
	- Updated figures produced by changepoint analysis

# sounDMR 2.2.2
	- Updated generate_methylframe() documentation

# sounDMR 2.2.1
	- Updated split_by_chromosome() to ignore scaffolds in the bed file
	- Bugfix for megaframe - was only saving unfiltered megaframe. Will now write out Megaframe after the NAs filter.

# sounDMR 2.2.0
	- find_DMR() local variable bugfix
	- Parameter fix for generate_methylframe() 
	- Minor updates to upstream functions
	- Added split_by_chunk that enables creating manageable chunks of the bedfile
	- updated split_by_chrom0some

# sounDMR 2.1.3
	- Bugfix for exporting max_read_depth parameter
	- minor enhancements to run generate_megaframe faster

# sounDMR 2.1.2

	- Minor readme updates
	- Minor updates to the functions
	- Bugfix in exporting split_by_chromosome function 

# sounDMR 2.1.1

* Package updated 
	-Includes bugfix for saving files with the prefix provided by the user
	-Includes bugfix for filtering the bedfile containing > max_read_depth paarmeter in the generate_methylframe the value for which can be provided by the user. 
	
# sounDMR 2.1.0

* Package updated to handle large bed files incase of whole genome sequencing
	-Includes fix for memory leaks
	-A built-in method to split large bed files by chromosome

# sounDMR 2.0.0

* Initial public release!
