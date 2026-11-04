# sounDMR <img src="man/figures/logo.png" align="right" height="120" alt="" />
For more information on the package, as well as examples of the types of downstream analyses that can be run using these results see Colicchio et al., 2023, [Scientific Reports](https://www.nature.com/articles/s41598-023-38356-7)

The following document contains a walkthrough of how to get ONT data into the appropriate format to run through sounDMR and a simple tutorial of the functions used to run a differential methylation analysis using sounDMR.

## What's new in 4.0.0
- **Downstream pipeline**: `run_dmr_downstream()` combines the significant-DMR tables for any design (one target gene for all treatments, one treatment with several genes, or different genes per treatment). It adds per-treatment methylation and read depth for every DMR, and writes summary tables, a BED file and figures.
- **`dmr_score4`**: rewards larger regions but gates out regions without a real signal; it is now the default in the example workflows.
- **Revised `boot_score()`**: tests the best DMR per context, the three combined, and the total of all DMRs near the target. It reports an excess score over the background, and a methylation effect (cytosine-equivalents) for comparing treatments between experiments.
- **Off-target treatments**: compare each DMR against treatments that targeted a different gene, as well as against the control.
- **Example workflows** in `inst/examples`: one script per experimental design, a step-through version and a per-individual analysis (see `inst/examples/README.md`).

See `NEWS.md` for details.

Note: Do not use the contents of the Docker. This is for Whole genome methylation analysis and is still in its beta stage. We are currently building a robust code to handle large data and will be releasing the full version soon.  


## Pre-requisites
We need to run Megalodon/Deep-Signal Plant/Bonito/Dorado to calculate methylation levels at each cytosine for each individual. This package technically works with any ONT bed file with a format as mentioned in [here](https://www.encodeproject.org/data-standards/wgbs/), in the **Description of bedMethyl file section**. If focusing on a small part fo the genome, these bed files must be subset using `bedtools intersect` to focus on specific regions of the genome.

```
bedtools intersect -a [ONT_methyl.bed] -b [target_regions.bed] -wa > [ONT_methyl_subset.bed]
```
More information on bedtools intersect can be found [here](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html).

In addition, these files must constain the cytosine contexts in the last column. This can be achieved using [modbam2bed](https://github.com/epi2me-labs/modbam2bed) for older hardware to split the bed file based on cytosine context. A simple `awk` command can then be used to add the cytosine context into an additional column before recombining into a single file. Alternatively this can be achieved with [modkit](https://github.com/nanoporetech/modkit) with the following line to convert into the proper format:
`cat ${unconverted_bed} | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8 ,$9, $10, $11, $4 }' |  sed 's/m,//g' |  sed 's/,0//g' > ${output_bed}`

## Installation
Make sure you have access to the latest version of the package and go through the installation steps
below. This can be done either by downloading the `tar.gz` file from
[github](https://github.com/SoundAg/sounDMR) or by cloning the project and moving to your local
environment or by installing directly from github via devtools. This only needs to be done once.

``` r
# Install directly from github
install.packages("devtools")
library(devtools)
install_github("SoundAg/sounDMR")

# Or from a downloaded tar.gz
install.packages("sounDMR_4.0.0.tar.gz", repos = NULL, type = "source")
```

## Usage
`library(sounDMR)` will load the core packages:

- [tidyverse](https://tidyverse.org)
- [reshape2](https://www.rdocumentation.org/packages/reshape2/versions/1.4.4)
- [changepoint](https://cran.r-project.org/web/packages/changepoint/changepoint.pdf)
- [lme4](https://cran.r-project.org/web/packages/lme4/lme4.pdf)
- [Formula](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/formula)
- [glmmTMB](https://cran.r-project.org/web/packages/glmmTMB/glmmTMB.pdf)

``` r
library(sounDMR)

#>Loading required package: changepoint
#>Loading required package: zoo
#>
#>Attaching package: ‘zoo’
#>
#>The following objects are masked from ‘package:base’:
#>
#>    as.Date, as.Date.numeric
#>
#>Successfully loaded changepoint package version 2.2.4
#> See NEWS for details of changes.
#>Loading required package: Formula
#>Loading required package: lme4
#>Loading required package: Matrix
#>Loading required package: reshape2
#>Loading required package: tidyverse
#>── Attaching packages ────────────────────────────────────────── tidyverse 1.3.2 ──
#>✔ ggplot2 3.4.0      ✔ purrr   0.3.5
#>✔ tibble  3.1.8      ✔ dplyr   1.0.10
#>✔ tidyr   1.2.1      ✔ stringr 1.4.1
#>✔ readr   2.1.3      ✔ forcats 0.5.2
#>── Conflicts ───────────────────────────────────────────── tidyverse_conflicts() ──
#>✖ tidyr::expand() masks Matrix::expand()
#>✖ dplyr::filter() masks stats::filter()
#>✖ dplyr::lag()    masks stats::lag()
#>✖ tidyr::pack()   masks Matrix::pack()
#>✖ tidyr::unpack() masks Matrix::unpack()
#>Loading required package: glmmTMB

```

## Data availability
Sample data is available on the [github](https://github.com/SoundAg/sounDMR) page which includes methyl bed files for 2 treated and control individuals each along with Experimental design file (file containing the experimental matrix data) and Gene coordinates file for the ZoomFrame.


## Workflow
This package follows the general workflow found in the example scripts in *inst/examples*. Start with *DMR_workflow.R*; *inst/examples/README.md* explains which script fits which design.

| Script | Use when |
|---|---|
| `DMR_workflow.R` | One target gene, shared by all treatments (default) |
| `DMR_workflow_multi_gene.R` | One treatment targets several genes |
| `DMR_workflow_treatment_specific_genes.R` | Different treatments target different genes |
| `DMR_workflow_stepwise.R` | Any design, run one treatment at a time with treatment-specific penalties, windows and cutoffs |
| `DMR_workflow_individual.R` | Per-plant analysis (each treated plant vs. all controls) around the target |

The steps are outlined below:

#### Create Methylframe
The first step is to generate Methylframe. This function calls 2 functions to generate megaframe and generate zoomframe depending on the type of analysis. Examples for both are provided in the Readme. Megaframe is a dataframe that includes combined methyl.bed information from all the samples in the analysis.Zoomframe is a version of the megaframe that includes additional columns like,
gene name or gene ID : There is flexibility to choose what column should be included in the dataframe using the Gene_column parameter, make sure to provide the exact same name;
Zeroth_pos : Position adjusted to have position 0 as ATG for the speicifc gene;
Gene : This can inclkude either Gene Id or Gene names for tracking during DMR analysis;
Zoom_co : This keeps track of whether the given position is within the gene or is an adaptive region.

```
# Find the Bedmethyl files within the working directory
methyl_bed <- list.files(path=".",pattern="*.bed")

```
##### Whole genome bedmethyl processing
For whole genome bed files, additional steps are necessary to ensure data is in the right format and that it is subset by chromsomes
to make sure that R doesn't run out of memory. 

```
#Split the large bed file by chromosomes
bedlist <- c()
for (i in 1:length(methyl_bed)) {
  beds <- split_by_chromosome(methyl_bed[i])
  bedlist[(length(bedlist) + 1)] <- list(beds)
}

#get the list of chromosome names to automatically used to grep the respective bedfiles
#or you can manually add a list of chromosomes chr_list <- c()
chrs_list <- unique(str_extract(bedlist,"ch0"))


#extract bedfile names based on the chromosome list.
#Note: For whole genome bed files, You can choose 1 at a time to avoid running out of memory
#update chrs_list[1] for the chromosome of interest
methyl_bed <- bedlist[grep(chrs_list[1],bedlist)]

```

Zoomframe is only created if and when the Geneco file is provided. Geneco file should include the following columns at the bare minimum :
Gene_Name | Chromosome | Low | High | Gene_length | Strand | Adapt_Low | Adapt_High .
Low and High refer to gene body coordinates in the 5'to 3' direction and Adapt_Low and Adapt_High refers to the coordinates for adaptive region around the gene. The previous step can be skipped in this case. 

```
# Import the Geneco file
Geneco <- read.table(file.choose(), header=TRUE, sep=",")


# create methylframe
#without gene info
Methylframe <- generate_methylframe(methyl_bed_list=methyl_bed, Sample_count = 0,
                                  Methyl_call_type="Dorado", filter_NAs = 0,
                                  gene_info = FALSE, gene_coordinate_file = NA, Gene_column=NA,
                                  target_info=FALSE,
                                  File_prefix="Sample")


#with gene info
Methylframe <- generate_methylframe(methyl_bed_list=All_methyl_beds, Sample_count = 0,
                                  Methyl_call_type="Dorado", filter_NAs = 0,
                                  gene_info = TRUE, gene_coordinate_file = Geneco, Gene_column='Gene_Name',
                                  target_info=TRUE,
                                  File_prefix="Sample")
```
Note: This function also saves an experimental design starter which would require you to include all the details before the analysis, along with megaframe and zoomframe(if gene_info=TRUE) in the current working directory.

#### Clean and Rearrange Data
This step cleans the data read into the environment and creates "long" formats
of the count of methylated and unmethylated reads as well as the percent
methylation for each individual.

```
experimental_design_df <- read.table(file.choose(), header=TRUE, sep=",")
dmr_obj <- create_dmr_obj(Methylframe, experimental_design_df)
```

#### Creating Methylation Summary Data
This section is involved in manipulating the methylation data to summarize
information based on specific variables such as the Percent methylation of each
individual for each cytosine. It also creates three columns for each individual:
the Z score, the change in methylation between the treatment and the control, and
the read depth for an individual. There is an optional `additional_summary_cols`
parameter that will allow us to include additional columns given a summary
statistic function such as `mean`, `sd`, or `var` and the name of the column on
which to run the summary statistics. Multiple tuples can be added into this
list to create additional summary columns. Note: the function name must be
a string.

```
methyl_summary <- create_methyl_summary(dmr_obj, control = 'C', treated = 'T',
                                        additional_summary_cols = list(c('sd', 'Group')))

# Option to subset methyl_summary
individuals_of_interest = unique(dmr_obj$experimental_design_df$Individual)
methyl_summary <- subset_methyl_summary(methyl_summary,
                                        individuals_to_keep = individuals_of_interest)
```

#### DMR Analysis
This section is where the model is run for our DMR analysis. We recommend using a
binomial model when comparing groups (e.g. control vs treated) and a beta-binomial
model when comparing individuals vs the control group. The user has the option to change the fixed
and random effects passed into the model depending on the goals of the experiment.

```
# Group Model
methyl_summary <- find_DMR(methyl_summary, dmr_obj, fixed = c('Group'),
                           random = c('Individual'), reads_threshold = 3,
                           control = 'C', model = 'binomial',
                           analysis_type = 'group')

# Individual Model
methyl_summary <- find_DMR(methyl_summary, dmr_obj, fixed = c('Group'),
                           random = c('Individual'), reads_threshold = 5,
                           control = 'C', model = 'beta-binomial',
                           analysis_type = 'individual')
```

#### Changepoint Analysis
The output from the previous step gets passed into our changepoint analysis
section. Here the user picks a column upon which to run the analysis. We
recommend using one of the Z_score columns created from the fixed effects provided
to the model.

```
# The genes of interest for the experiment
target_genes <- unique(dmr_obj$ZoomFrame_filtered$Gene)

methyl_summary_cg <- changepoint_analysis(methyl_summary, CG_penalty = 9,
                                       CHG_penalty = 4, CHH_penalty = 7,
                                       target_genes = target_genes,
                                       save_plots = F,
                                       z_col = "Z_GroupT_small")
```
With `save_plots = TRUE`, `plot_prefix` (e.g. the treatment name) and `plot_dir` keep plots from different treatments apart, e.g. `plot_prefix = "T2", plot_dir = "changepoint_plots"`.

#### Sound score
This function takes the changepoint regions from the previous step, built on a specific test statistic, and creates an aggregated region summary with statistics for each region. A "sound score" is then computed for each region to measure the strength of a differentially methylated region (DMR). The scores are:

```
dmr_score  = Count^(1/3) * (|mean Z| * |mean change|)^(1/2)
dmr_score2 = Count^(1/3) * (|mean Z| * |asin(sqrt(change/100 + control/100)) - asin(sqrt(control/100))|)^(1/2)
dmr_score3 = Count^(1/3) * |mean Z|^(1/2) * |mean change|^(1/2) * |max Z|^(1/2)
dmr_score4 = scale * Count^e * |mean Z|^(1/2) * |mean change|^(1/2) * |max Z|^(1/2)
```

`dmr_score4` is recommended and is the default in the example workflows.
- The count exponent `e` (`score4_count_exponent`, default 0.5) rewards larger regions more than `dmr_score3` does.
- A signal gate sets the score to 0 for regions with |mean Z| < `score4_min_mean_z` (0.175) or |mean change| < `score4_min_mean_dif` (2 points). This stops large regions with no real signal from scoring highly.
- `score4_scale` (0.7) keeps it on about the same scale as `dmr_score3`.
- Whether a region passed the gate is in `dmr_score4_signal`.

The user can also give a custom DMR score function, e.g. `UserFunction = (Ag_Groups$Count)^(1/3)`. It can use any column of the region summary, referred to as `Ag_Groups`. This function also produces a plot of the distribution of DMR scores relative to their percentile.

```
DMR_score <- sound_score(changepoint_OF = methyl_summary_cg,
                         Statistic = "Z_GroupT_small",
                         Per_Change = "Treat_V_Control",
                         other_columns = c("Control"),
                         score4_count_exponent = 0.5, score4_scale = 0.7,
                         score4_min_mean_z = 0.175, score4_min_mean_dif = 2)
```

#### Boot score
This function tests an a priori hypothesis: that methylation changed around a target region, such as an oligo target window. It compares the DMRs near that region with random windows elsewhere.

Each region of the target gene gets a distance weight. The weight is 1 if the region overlaps the target window and falls by 0.00006 per bp of distance, reaching 0 beyond 15 kb. A region's adjusted score is its weight × its region score. The same is done for `nboots` random windows of the same width, placed in the other genes. The output reports:

| Test | What is compared | Output |
|---|---|---|
| Best DMR, per context | Best adjusted score near the target vs. the best in each background window | p-value, z, evidence |
| Best DMRs combined | Fisher combination of the three context p-values, calibrated on the same background windows | p-value |
| Total of all DMRs | Summed adjusted score of all regions above `overall_cutoff`, all contexts | total, background mean / 95th percentile, `Excess_Score` (total − background mean), z, p-value |
| Methylation effect | Sum of weight × cytosines × mean change / 100 for the same regions, in cytosine-equivalents (10 cytosines each +50 points = 5) | total, background mean / 95th percentile, `Excess_Meth_Effect` |

p-values are the share of background windows scoring at least as high. Evidence labels are strong (p < 0.01), significant (p < 0.05) or suggestive (p < 0.1). `direction_DMR = "negative"` scores loss of methylation. `exclude_genes` drops other target genes of the same treatment from the background. The function prints a summary table and returns it as `$summary`, together with `$plot` (background distributions with the target marked) and `$region_plot` (adjusted scores along the target gene).

```
DMR_boot <- boot_score(sound_score_obj = DMR_score, target_gene = "AT1G01640",
                       target_start = -150, target_end = 3000,
                       scoring_col_name = "dmr_score4", overall_cutoff = 8)
DMR_boot$summary

# Recompute the total-of-DMRs test at another cutoff
DMR_boot <- boot_overall(DMR_boot, cutoff = 10, verbose = TRUE)
```

#### Significant DMRs
`get_sig_dmrs()` returns the target-gene regions whose adjusted score is above a cutoff. Use the same cutoff as `overall_cutoff`, so the bootstrap total matches the summed score of the DMRs you keep.

```
DMR_boot_T2_sig <- get_sig_dmrs(DMR_boot, min_score = 8)
```

#### Downstream summaries
`run_dmr_downstream()` combines the significant-DMR tables of all treatments. The treatment comes from each table's name (`DMR_boot_<TX>_sig`) and the gene from its `Gene` column, so any design works. It then adds per-treatment methylation and read depth for every DMR and writes everything to `out_dir`.

```
res <- run_dmr_downstream(sig_list = mget(ls(pattern = "^DMR_boot_.*_sig$")),
                          dmr_obj, control = "Water",
                          off_targets = list(T2 = c("T6", "T7")),   # optional
                          boots = dmr_boots,                        # optional: bootstrap results
                          out_dir = "db210_downstream", prefix = "db210")
```

| File | Contents |
|---|---|
| `_dmrs.csv` | One row per DMR: scores, methylation effect, target / control (/ off-target) % methylation and read depth |
| `_summary.csv`, `_summary_by_gene.csv` | Per treatment (× gene): summed score, DMR counts, maximum score and methylation effect by context; with `boots`, the bootstrap p-values, z, backgrounds and excess scores |
| `_bootstrap.csv` | Bootstrap results per treatment × gene |
| `_dmr_group_stats.csv` | One row per DMR × treatment group |
| `_sig_table_mapping.csv` | Which sig table → treatment → gene(s) |
| `_dmrs.bed` | BED4 of the significant DMRs |
| `_score_summary.pdf` | Summed score per treatment by context, with background mean / 95th-percentile ticks, excess and p |
| `_dmr_map.pdf`, `_dmr_methylation.pdf` | DMR positions, and methylation per treatment group for each DMR |

**Off-target treatments** are treatments that targeted a different gene. Set them with `off_targets`. Each DMR then also gets its difference from the off-target treatments, an `Off_Target_Fraction` and an `Off_Target_Adjusted_Score`. This shows whether a change is specific to the treatments that targeted the gene.

**Comparing treatments.** Within an experiment, use the scores and p-values (`Excess Total Score`, `P Total Score`). Between experiments, use `Excess Meth Effect`. Scores are built from Z statistics, so they grow with read depth and number of plants, and they depend on the changepoint penalties and the scoring function. The methylation effect is in units of methylation change and is not affected by either. Keep the target window and cutoff the same between experiments. `inst/examples/README.md` lists every summary column.
