methyl_bed <- list.files(path=".",pattern="*.bed")

Geneco <- read.table(file.choose(), header=TRUE)

Methylframe <- generate_methylframe(methyl_bed_list=methyl_bed, Sample_count = 0,
                                    Methyl_call_type="Modkit", filter_NAs = 6,
                                    gene_info = TRUE, gene_coordinate_file = Geneco, Gene_column='Gene_Name',
                                    target_info=FALSE,
                                    File_prefix="Sample")

experimental_design_df <- read.table(file.choose(), header=TRUE, sep=",")

head(experimental_design_df)

dmr_obj <- create_dmr_obj(Methylframe, experimental_design_df)

methyl_summary <- create_methyl_summary(dmr_obj, control = 'C', treated = 'T',
                                        additional_summary_cols = list(c('mean', 'Treatment')))


methyl_summary <- find_DMR(methyl_summary, dmr_obj, fixed = c('Treatment'),
random = c('Individual'), reads_threshold = 3,
control = 'Water', model = 'binomial',
analysis_type = 'group', target_variable = "Treatment")

target_genes <- unique(dmr_obj$ZoomFrame_filtered$Gene)


methyl_summary$T6_mean_dif <- methyl_summary$T6_mean- methyl_summary$Control


#For use in case where there are off target groups
#methyl_summary$T6_mean-((methyl_summary$T6_mean+methyl_summary$T6_mean+methyl_summary$T6_mean+methyl_summary$T6_mean)/4)->methyl_summary$T6_off

methyl_summary <- changepoint_analysis(methyl_summary, CG_penalty = 3,
                                            CHG_penalty =3, CHH_penalty =5,
                                            target_genes = target_genes,
                                            save_plots = TRUE,
                                            z_col = "Z_TreatmentT6_small")


#methyl_summary_tiny <- changepoint_analysis(methyl_summary_tiny, CG_penalty = 1200,
                                            #CHG_penalty =600, CHH_penalty =2400,
                                            #target_genes = target_genes,
                                            #save_plots = TRUE,
                                            #z_col = "AGL24_Dif")
#In case of multple different targets within one experiment carry over the "_off" columns here

DMR_score_T6 <- sound_score(changepoint_OF = methyl_summary ,
                             Statistic = "Z_TreatmentT6_small",
                             Per_Change = "T6_mean_dif", CF = F,
                             other_columns=c("Treat_V_Control"),
                             UserFunction = NA)

#add Statistic to keep col, update summarise at, 

#DMR_score_AGL24_Delta <- sound_score(changepoint_OF = methyl_summary_tiny ,
                                   #Statistic = "AGL24_Dif",
                                   #Per_Change = "AGL24_Dif", CF = F,
                                   #other_columns=c("H20_mean"),
                                   #UserFunction = NA)


DMR_boot_T6_PAL1 <- boot_score(sound_score_obj = DMR_score_T6,
                           target_gene = "LsPAL1", scoring_col_name="dmr_score3", target_start = -2000, target_end = 1000 )

DMR_boot_T6_PAL2 <- boot_score(sound_score_obj = DMR_score_T6,
                               target_gene = "LsPAL2", scoring_col_name="dmr_score3", target_start = -2000, target_end = 1000 )


#DMR_boot_AGL24_Delta <- boot_score(sound_score_obj = DMR_score_AGL24_Delta,
                                 #target_gene = "AGL24", scoring_col_name="dmr_score_noZ", target_start = -1382, target_end = -8515 )

DMR_boot_T6_PAL1$target_rs[DMR_boot_T6_PAL1$target_rs$adjusted_soundscore>4,]->DMR_boot_T6_PAL1_sig

treatment_plots(sig_dmrs = DMR_boot_T6_PAL1_sig , dmr_score_file = DMR_score_T6 , prefix = "T6_PAL1_" )

DMR_boot_T6_PAL2$target_rs[DMR_boot_T6_PAL2$target_rs$adjusted_soundscore>4,]->DMR_boot_T6_PAL2_sig

treatment_plots(sig_dmrs = DMR_boot_T6_PAL2_sig , dmr_score_file = DMR_score_T6 , prefix = "T6_PAL2_" )



# Get names of data frames ending with "_data"
data_frames_to_combine_names <- ls(pattern = "_sig$")

# Retrieve the data frame objects
data_frames_list <- mget(data_frames_to_combine_names)

new_colnames <- colnames(data.frame(mget(data_frames_to_combine_names[1])))

list_of_dfs_renamed <- lapply(data_frames_list, function(df) {
  colnames(df) <- new_colnames
  return(df)
})


# Combine them using do.call and rbind
combined_df <- rbindlist(list_of_dfs_renamed, idcol = "source_df")


combined_df[,c(3,6,7,1,9, 31)]->cdf_bed
sub("DMR_boot_", "",cdf_bed$source_df)->cdf_bed$source_df
sub("_sig", "",cdf_bed$source_df)->cdf_bed$source_df
round(cdf_bed$DMR_boot_T2_PAL1_sig.adjusted_soundscore, digits = 2)->cdf_bed$DMR_boot_T2_PAL1_sig.adjusted_soundscore
paste(cdf_bed$source_df, cdf_bed$DMR_boot_T2_PAL1_sig.CX, cdf_bed$DMR_boot_T2_PAL1_sig.adjusted_soundscore, sep = "_")->cdf_bed$dmr_title
cdf_bed[,c(1,2,3,7)]->cdf_bed
write.table(cdf_bed, file="db172_dmr.bed")

#Rename Summary Table Columns
names(combined_df) <- gsub("DMR_boot_T2_PAL1_sig.", "", names(combined_df))
names(combined_df) <- gsub("_TreatmentT2_small", "", names(combined_df))
#names(combined_df) <- gsub("T2_mean_dif", "Dif_Vs_OffTarget", names(combined_df))
gsub("DMR_boot_", "",combined_df$source_df)->combined_df$source_df
gsub("_.*", "",combined_df$source_df)->combined_df$source_df

DMR_Info(combined_DMRs = combined_df, dmr_obj_file=dmr_obj)->combined_df


library(modelsummary)

summary_table<-datasummary(factor(Gene)*factor(source_df)*(adjusted_soundscore) ~ (sum + N + Max)*factor(CX),
data = combined_df, output = ".csv")

size_table<-datasummary(factor(Gene)*factor(source_df)*(Count) ~ (sum + N + Max),
                           data = combined_df, output = ".csv")


cbind(summary_table[,-3],size_table$sum)->summary_table
summary_table$Total_Score<- (as.numeric(summary_table[,3])+as.numeric(summary_table[,4])+as.numeric(summary_table[,5]))
summary_table$Description<-NA
summary_table$Oligo_Duplexes<-NA
summary_table<-summary_table[,c(1,2,14,15,13,3:12)]
colnames(summary_table)<-c("Gene", "Treatment", "Description","Oligo Duplexes", "Total_Score", "CG Score", "CHG Score", "CHH Score", "# CG DMRs", "# CHG DMRs", "# CHH DMRs", "Max CG Score", "Max CHG Score", "Max CHH Score", "Cytosines Covered")


write.csv(summary_table, "db_172_summary.csv", row.names = FALSE)
write.csv(combined_df, file="db_172_DMRs.csv")


