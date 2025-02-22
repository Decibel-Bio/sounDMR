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



library(sounDMR)

#-----------------Part 1 : ONT data clean up and standardization----------------

#-------------------------------
# Read In Gene Co-ordinates File
#-------------------------------
# If available, it is necessary for the gene coordinates file to be in the below 
# format to ensure the code works.  
# Chromosome|Low|High|Gene_name|Strand|Gene_length|Adapt_Low|Adapt_High

Geneco <- read.table(file.choose(), header=TRUE, sep=",")

#-------------------------
# Read In Methyl.bed Files
#-------------------------
# The below command works only if all the bed files that you want to work with 
# are in the working directory.
# Make sure to change the pattern based on your methyl bed file
methyl_bed <- list.files(path=".",pattern="*.bed")

#-------------------------------------
# Whole Genome Bed File Processing Only
#-------------------------------------

# NOTE: Only applicable for whole genome bed files; jump to the next section if 
# you have used adaptive sequencing
# Split the large bed file by chromosomes

# Get the list of chromosome names to automatically used to grep the respective 
# bedfiles or you can manually add a list of chromosomes chr_list <- c()

# Extract bedfile names based on the chromosome list.
# Note: For whole genome bed files, You can choose 1 at a time to avoid running 
# out of memory
# Note: update chrs_list[1] for the chromosome of interest

#-------------------------
# Create Methylframe 
#-------------------------
# If gene_info is false in the below parameter then this returns a megaframe and 
# if True it returns the zoomframe
Methylframe <- generate_methylframe(methyl_bed_list=methyl_bed, 
                                    Sample_count = 0, 
                                    Methyl_call_type="Dorado", 
                                    filter_NAs = 1, 
                                    max_read_depth=100,
                                    gene_info = FALSE, 
                                    gene_coordinate_file = NA, 
                                    Gene_column=NA,
                                    target_info=FALSE, 
                                    File_prefix="Sample")



# Note: The above function creates The experimental_design starter doesn't have 
# any information with respect to treatments, rounds etc. 
# Make sure to add it for DMR analysis

#----------------------------Part 2 : DMR analysis------------------------------

#------------
# Clean Data
#------------
# Note: Update the created "Sample_Experimental_design_starter.csv" file to
# include the other information
dmr_obj <- create_dmr_obj(Methylframe, experimental_design_df)

#-----------------------------------------------
# Creating Differential Methylation Output File
#-----------------------------------------------

methyl_summary <- create_methyl_summary(dmr_obj, 
                                        control = 'C', 
                                        treated = 'T',
                                        additional_summary_cols = list(
                                          c('sd', 'Group')
                                          ))

# Option to subset methyl_summary
# Note: code is currently set up to include all individuals, change as needed
individuals_of_interest = unique(dmr_obj$experimental_design_df$Individual)
methyl_summary <- subset_methyl_summary(methyl_summary, 
                                        individuals_to_keep=individuals_of_interest)

#--------------------
# Group DMR Analysis
#--------------------

# Run the Group DMR analysis
methyl_summary <- find_DMR(methyl_summary, dmr_obj, fixed = c('Group'), 
                           random = c('Individual'), reads_threshold = 3, 
                           control = 'C', model = 'binomial', 
                           analysis_type = 'group')

#----------------
# Individual DMR
#----------------
# Run the Individual DMR analysis
methyl_summary <- find_DMR(methyl_summary, dmr_obj, fixed = c('Group'), 
                           random = c('Individual'), reads_threshold = 5, 
                           control = 'C', model = 'beta-binomial', 
                           analysis_type = 'individual')

#----------------------
# Changepoint Analysis
#----------------------
# Get the potential column names to run changepoint analysis on
changepoint_cols = find_changepoint_col_options(methyl_summary)

# The target genes of interest
target_genes <- unique(dmr_obj$ZoomFrame_filtered$Gene)

# Run the changepoint_analysis function
# Note: when whole_genome = FALSE, the target_genes need to be a non-empty list
# Note: changepoint analysis requires methyl_summary to be sorted by Chromosome
# and position
methyl_summary <- changepoint_analysis(methyl_summary, 
                                       CG_penalty = 5, 
                                       CHG_penalty = 7, 
                                       CHH_penalty = 6, 
                                       target_genes = c(),
                                       save_plots = TRUE,
                                       z_col = "Z_GroupT_small", 
                                       whole_genome = TRUE)

#----------------------
# DMR score rendering
#----------------------

DMR_score <- sound_score(changepoint_OF = methyl_summary, 
                         Statistic = changepoint_cols[1], 
                         Per_Change = "Treat_V_Control", CF = T,
                         other_columns=c("Control", "Estimate_GroupT_small"),
                         UserFunction = NA)

# Only run bootscore if gene info is available
DMR_boot_score <- boot_score(sound_score_obj = DMR_score, 
                             target_gene = "AT1G01640", 
                             scoring_col_name="dmr_score2")
