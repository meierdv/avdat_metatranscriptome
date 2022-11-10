library("DESeq2")
require(data.table)

All_Feature_counts <- read.csv("Avdat_featureCounts.tsv", 
                               sep="\t", 
                               header=T, 
                               quote="",
                               check.names=F)
All_Feature_counts <- All_Feature_counts[,c(1,2,6,7:ncol(All_Feature_counts))]

#Using the opportunity to rename samples
colnames(All_Feature_counts) <- c("Geneid","Contig","Length","C_0","E_0","F_0","C_0.25","E_0.25","F_0.25","C_0.5","E_0.5","F_0.5","C_3","E_3","F_3","C_6","E_6","F_6","C_12","E_12","F_12","C_29","E_39","F_39","C_55","E_55","F_55")

#Fragments per kbp
FPK <- All_Feature_counts[,3:ncol(All_Feature_counts)] / (All_Feature_counts[,2]/1000)

#Transcripts per kbp per million
TPM <- sweep(FPK,2, colSums(FPK)/1000000, "/")

TPM$Geneid <- All_Feature_counts$Geneid

Annotation <- read.csv("Avdat_crust_all_annotations.tsv", 
                       header=T, 
                       sep="\t", 
                       quote="",
                       check.names=F)

write.table(merge(Annotation, TPM, by="Geneid", all.y=T), 
            file="Avdat_featureCounts_TPM.tsv", 
            sep="\t", 
            row.names=F, 
            quote=F)

#generate subtotals for columns 5 - 32. function to be applied is "sum"
Transcripts_per_bin <- cube(merge(Annotation[,c("Bin","Geneid")], TPM, by="Geneid", all.y=T), 
                            lapply(.SD, sum), 
                            by="Bin", 
                            .SDcol=c(3:26))

#generate subtotals for columns 5 - 32. function to be applied is "count non-zeroes"
Transcribed_genes_per_bin <- cube(merge(Annotation[,c("Bin","Geneid")], TPM, by="Geneid", all.y=T), 
                                        lapply(.SD, function(c)sum(c!=0)), 
                                  by="Bin", 
                                  .SDcol=c(3:26))

write.table(Transcripts_per_bin, file="Avdat_TPMs_per_bin.tsv", sep="\t", row.names=F, quote=F)
write.table(Transcribed_genes_per_bin, file="Avdat_transcribed_genes_per_bin.tsv", sep="\t", row.names=F, quote=F)


#---------------------------------------------------------------------------------------------------
#Generating TPM values normalized per bin
#----------------------------------------------------------------------------------------------------

system("bash Subset_featureCounts_table.sh", wait=T)

for(i in read.table("Bin_list.txt")[,1]) {
  
  #Importing the feature Counts
  countsfile <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts.tsv", sep=""))
  FeatureCounts <- read.table(countsfile, 
                              sep="\t", 
                              header=F, 
                              quote="")
 
  colnames(FeatureCounts) <- colnames(All_Feature_counts)
  FPK <- FeatureCounts[,4:27] / (FeatureCounts[,3]/1000)
  TPM <- sweep(FPK,2, colSums(FPK)/1000000, "/")
  TPM$Geneid <- FeatureCounts$Geneid
  TPM$Contig <- FeatureCounts$Contig
  TPM <- TPM[,c(25,26,1:24)]

  output <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts_TPMperbin.tsv", sep=""))
  
  write.table(
    merge(Annotation, TPM, by="Geneid", all.y=T),
    file=output, 
    sep="\t", 
    row.names=F, 
    quote=F)
}


#---------------------------------------------------------------------------------------------------------
#Differential expression analysis
#---------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------------
#Analyzing differential expression between e.g. 0 min and 15 min timepoints of the rehydration experiment
#---------------------------------------------------------------------------------------------------------------------------------
#DeSeq2 compares one category vs another. Therefore a subsetting of the Metadata and Featurecounts tables is necessary.
#Importing metadata (subsetting to compare 0 min and 15 min timepoints)
rehydration_metadata_sub <- read.table("Metadata.csv", 
                                     sep="\t", 
                                     header=T, 
                                     row.names=1, 
                                     quote="")[c(1:6),]
#Bin_list_deseq.txt is the list of bins fulfilling the criteria for Deseq2 analysis. Otherwise, Deseq2 will abort.

for(i in read.table("Bin_list_deseq.txt")[,1]) {
  
  #Importing the feature Counts
  countsfile <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts.tsv", sep=""))
  
  FeatureCounts <- read.table(countsfile, 
                              row.names=1, 
                              sep="\t", 
                              header=F, 
                              quote="")[,c(4:9)]
  
  colnames(FeatureCounts) <- c("C_0","E_0","F_0","C_0.25","E_0.25","F_0.25")
  
  #The TPM table for the MAG will later be merged with the Deseq results, so that transcription profiles of significantly diff. expressed genes can be checked on the spot
  TPM_file <- file.path(paste("Transcripts_per_bin/", i, "_featureCounts_TPMperbin.tsv", sep=""))
  TPM_table <- read.csv(TPM_file, 
                        sep="\t", 
                        header=T, 
                        quote="")
  
  #Subsetting the count table in order to compare two conditions against each other. 
  #In this case 0 min and 15 min after rehydration
  
  FeatureCounts_sub <- FeatureCounts[,1:6]
  
  dds <- DESeqDataSetFromMatrix(countData = FeatureCounts_sub,
                                colData = rehydration_metadata_sub,
                                design= ~ Timepoint)
  
  dds <- DESeq(dds)
  #resultsNames(dds) # lists the coefficients
  res <- as.data.frame(results(dds, name="Timepoint_T00.25_vs_T00"))
  setDT(res, keep.rownames = "Geneid")
  res <- transform(res, Geneid = as.numeric(Geneid))
  
  Deseq_results <- merge(res, TPM_table, 
          by="Geneid", 
          all.x=T)
  
  output <- file.path(paste("./Diff_expression/", i, "_00.25_00.tsv", sep=""))
  
  write.table(Deseq_results, 
              file=output, 
              sep="\t",
              row.names=F,
              quote=F)
  
}


#---------------------------------------------------------------------------------------------------------------------------------
#Analyzing differential expression between all "dry" and all "hydrated" timepoints of the rehydration experiment
#---------------------------------------------------------------------------------------------------------------------------------

#Importing metadata (subsetting to compare all dry and all hydrated timepoints)
rehydration_metadata_sub <- read.table("Metadata.csv", 
                                     sep="\t", 
                                     header=T, 
                                     row.names=1, 
                                     quote="")[c(1:18,20:24),]

for(i in read.table("Bin_list_deseq.txt")[,1]) {
  
  #Importing the feature Counts
  countsfile <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts.tsv", sep=""))
  
  FeatureCounts <- read.table(countsfile, 
                              row.names=1, 
                              sep="\t", 
                              header=F, 
                                quote="")[,c(3:20,22:26)]
  
  colnames(FeatureCounts) <- c("C_0","E_0","F_0","C_0.25","E_0.25","F_0.25","C_0.5","E_0.5","F_0.5","C_3","E_3","F_3","C_6","E_6","F_6","C_12","E_12","F_12","E_39","F_39","C_55","E_55","F_55")
  
  
  TPM_file <- file.path(paste("Transcripts_per_bin/", i, "_featureCounts_TPMperbin.tsv", sep=""))
  TPM_table <- read.csv(TPM_file, 
                        sep="\t", 
                        header=T, 
                        quote="")
  
  
  #Subsetting the count table in order to compare two conditions against each other. 
  #In this case all hydrated vs all dry
  
  FeatureCounts_sub <- FeatureCounts[,1:23]
  
  dds <- DESeqDataSetFromMatrix(countData = FeatureCounts_sub,
                                colData = rehydration_metadata_sub,
                                design= ~ Phase1)
  
  dds <- DESeq(dds)
  #resultsNames(dds) # lists the coefficients
  res <- as.data.frame(results(dds, name="Phase1_hydrated_vs_dry"))
  setDT(res, keep.rownames = "Geneid")
  res <- transform(res, Geneid = as.numeric(Geneid))
  
  #Deseq_results <- merge(res, TPM_table, 
  #                       by="Geneid", 
  #                       all.x=T)
  
  output <- file.path(paste("./Diff_expression/", i, "_allhydrated_alldry.tsv", sep=""))
  
  write.table(res, 
              file=output, 
              sep="\t",
              row.names=F,
              quote=F)
  
}


#---------------------------------------------------------------------------------------------------------------------------------
#Analyzing differential expression between "dry" and "early hydrated" timepoints of the rehydration experiment
#---------------------------------------------------------------------------------------------------------------------------------

#Importing metadata (subsetting to compare eraly hydrated and dry timepoints)
rehydration_metadata_sub <- read.table("Metadata.csv", 
                                     sep="\t", 
                                     header=T, 
                                     row.names=1, 
                                     quote="")[c(1:9,20:24),]

for(i in read.table("Bin_list_deseq.txt")[,1]) {
  
  #Importing the feature Counts
  countsfile <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts.tsv", sep=""))
  
  FeatureCounts <- read.table(countsfile, 
                              row.names=1, 
                              sep="\t", 
                              header=F, 
                              quote="")[,c(3:20,22:26)]
  
  colnames(FeatureCounts) <- c("C_0","E_0","F_0","C_0.25","E_0.25","F_0.25","C_0.5","E_0.5","F_0.5","C_3","E_3","F_3","C_6","E_6","F_6","C_12","E_12","F_12","E_39","F_39","C_55","E_55","F_55")
  
  
  TPM_file <- file.path(paste("Transcripts_per_bin/", i, "_featureCounts_TPMperbin.tsv", sep=""))
  TPM_table <- read.csv(TPM_file, 
                        sep="\t", 
                        header=T, 
                        quote="")
  
  
  #Subsetting the count table in order to compare two conditions against each other. 
  #In this case dry and earlyhydrated
  
  FeatureCounts_sub <- FeatureCounts[,c(1:9,19:23)]
  
  
  dds <- DESeqDataSetFromMatrix(countData = FeatureCounts_sub,
                                colData = rehydration_metadata_sub,
                                design= ~ Phase2)
  
  dds <- DESeq(dds)
  #resultsNames(dds) # lists the coefficients
  res <- as.data.frame(results(dds, name="Phase2_earlyhydrated_vs_dry"))
  setDT(res, keep.rownames = "Geneid")
  res <- transform(res, Geneid = as.numeric(Geneid))
  
  #Deseq_results <- merge(res, TPM_table, 
  #                       by="Geneid", 
  #                       all.x=T)
  
  output <- file.path(paste("./Diff_expression/", i, "_earlyhydrated_dry.tsv", sep=""))
  
  write.table(res, 
              file=output, 
              sep="\t",
              row.names=F,
              quote=F)
  
}

#---------------------------------------------------------------------------------------------------------------------------------
#Analyzing differential expression between "early hydrated" and "hydrated" timepoints of the rehydration experiment
#---------------------------------------------------------------------------------------------------------------------------------

#Importing metadata (subsetting to compare early hydrated and hydrated timepoints)
rehydration_metadata_sub <- read.table("Metadata.csv", 
                                     sep="\t", 
                                     header=T, 
                                     row.names=1, 
                                     quote="")[c(4:9,10:18),]

for(i in read.table("Bin_list_deseq.txt")[,1]) {
  
  #Importing the feature Counts
  countsfile <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts.tsv", sep=""))
  
  FeatureCounts <- read.table(countsfile, 
                              row.names=1, 
                              sep="\t", 
                              header=F, 
                              quote="")[,c(3:20,22:26)]
  
  colnames(FeatureCounts) <- c("C_0","E_0","F_0","C_0.25","E_0.25","F_0.25","C_0.5","E_0.5","F_0.5","C_3","E_3","F_3","C_6","E_6","F_6","C_12","E_12","F_12","E_39","F_39","C_55","E_55","F_55")
  
  
  TPM_file <- file.path(paste("Transcripts_per_bin/", i, "_featureCounts_TPMperbin.tsv", sep=""))
  TPM_table <- read.csv(TPM_file, 
                        sep="\t", 
                        header=T, 
                        quote="")
  
  
  #Subsetting the count table in order to compare two conditions against each other. 
  #In this case hydrated and earlyhydrated
  
  FeatureCounts_sub <- FeatureCounts[,c(4:9,10:18)]
  
  
  dds <- DESeqDataSetFromMatrix(countData = FeatureCounts_sub,
                                colData = rehydration_metadata_sub,
                                design= ~ Phase2)
  
  dds <- DESeq(dds)
  #resultsNames(dds) # lists the coefficients
  res <- as.data.frame(results(dds, name="Phase2_hydrated_vs_earlyhydrated"))
  setDT(res, keep.rownames = "Geneid")
  res <- transform(res, Geneid = as.numeric(Geneid))
  
  #Deseq_results <- merge(res, TPM_table, 
  #                       by="Geneid", 
  #                       all.x=T)
  
  output <- file.path(paste("./Diff_expression/", i, "_hydrated_earlyhydrated.tsv", sep=""))
  
  write.table(res, 
              file=output, 
              sep="\t",
              row.names=F,
              quote=F)
  
}

#---------------------------------------------------------------------------------------------------------------------------------
#Analyzing differential expression between "hydrated" and "dry" timepoints of the rehydration experiment
#---------------------------------------------------------------------------------------------------------------------------------

#Importing metadata (subsetting to compare hydrated and dry timepoints)
rehydration_metadata_sub <- read.table("Metadata.csv", 
                                     sep="\t", 
                                     header=T, 
                                     row.names=1, 
                                     quote="")[c(1:3,10:18,20:24),]

for(i in read.table("Bin_list_deseq.txt")[,1]) {
  
  #Importing the feature Counts
  countsfile <- file.path(paste("./Transcripts_per_bin/", i, "_featureCounts.tsv", sep=""))
  
  FeatureCounts <- read.table(countsfile, 
                              row.names=1, 
                              sep="\t", 
                              header=F, 
                              quote="")[,c(3:20,22:26)]
  
  colnames(FeatureCounts) <- c("C_0","E_0","F_0","C_0.25","E_0.25","F_0.25","C_0.5","E_0.5","F_0.5","C_3","E_3","F_3","C_6","E_6","F_6","C_12","E_12","F_12","E_39","F_39","C_55","E_55","F_55")
  
  
  TPM_file <- file.path(paste("Transcripts_per_bin/", i, "_featureCounts_TPMperbin.tsv", sep=""))
  TPM_table <- read.csv(TPM_file, 
                        sep="\t", 
                        header=T, 
                        quote="")
  
  
  #Subsetting the count table in order to compare two conditions against each other. 
  #In this case hydrated and dry
  
  FeatureCounts_sub <- FeatureCounts[,c(1:3,10:23)]
  
  
  dds <- DESeqDataSetFromMatrix(countData = FeatureCounts_sub,
                                colData = rehydration_metadata_sub,
                                design= ~ Phase2)
  
  dds <- DESeq(dds)
  #resultsNames(dds) # lists the coefficients
  res <- as.data.frame(results(dds, name="Phase2_hydrated_vs_dry"))
  setDT(res, keep.rownames = "Geneid")
  res <- transform(res, Geneid = as.numeric(Geneid))
  
  #Deseq_results <- merge(res, TPM_table, 
  #                       by="Geneid", 
  #                       all.x=T)
  
  output <- file.path(paste("./Diff_expression/", i, "_hydrated_dry.tsv", sep=""))
  
  write.table(res, 
              file=output, 
              sep="\t",
              row.names=F,
              quote=F)
  
}

    

#---------------------------------------------------------------------------------------------------------------------------------
# Combining tables into one giant table
#---------------------------------------------------------------------------------------------------------------------------------
system("tail -n +2 -q Diff_expression/*_00.25_00.tsv > Diff_expression/Diff_expression_00.25_00.tsv", wait=T)
Log_pvalue_0_0.25 <- read.csv("Diff_expression/Diff_expression_00.25_00.tsv", 
                              header=T, 
                              sep="\t",
                              quote="",
                              check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_0_0.25) <- c("Geneid", "Log2chnage_00_0.25", "padj_00_0.25")

system("tail -n +2 -q Diff_expression/*_00.5_00.25.tsv > Diff_expression/Diff_expression_00.5_00.25.tsv", wait=T)
Log_pvalue_0.25_0.5 <- read.csv("Diff_expression/Diff_expression_00.5_00.25.tsv", 
                                header=T, 
                                sep="\t",
                                quote="",
                                check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_0.25_0.5) <- c("Geneid", "Log2chnage_0.25_0.5", "padj_0.25_0.5")

system("tail -n +2 -q Diff_expression/*_03_00.5.tsv > Diff_expression/Diff_expression_03_00.5.tsv", wait=T)

Log_pvalue_0.5_03 <- read.csv("Diff_expression/Diff_expression_03_00.5.tsv", 
                              header=T, 
                              sep="\t",
                              quote="",
                              check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_0.5_03) <- c("Geneid", "Log2chnage_0.5_03", "padj_0.5_03")

system("tail -n +2 -q Diff_expression/*_06_03.tsv > Diff_expression/Diff_expression_06_03.tsv", wait=T)
Log_pvalue_03_06 <- read.csv("Diff_expression/Diff_expression_06_03.tsv", 
                             header=T, 
                             sep="\t",
                             quote="",
                             check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_03_06) <- c("Geneid", "Log2chnage_03_06", "padj_03_06")

system("tail -n +2 -q Diff_expression/*_12_06.tsv > Diff_expression/Diff_expression_12_06.tsv", wait=T)
Log_pvalue_06_12 <- read.csv("Diff_expression/Diff_expression_12_06.tsv", 
                             header=T, 
                             sep="\t",
                             quote="",
                             check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_06_12) <- c("Geneid", "Log2chnage_06_12", "padj_06_12")

system("tail -n +2 -q Diff_expression/*_39_12.tsv > Diff_expression/Diff_expression_39_12.tsv", wait=T)
Log_pvalue_12_39 <- read.csv("Diff_expression/Diff_expression_39_12.tsv", 
                             header=T, 
                             sep="\t",
                             quote="",
                             check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_12_39) <- c("Geneid", "Log2chnage_12_39", "padj_12_39")

system("tail -n +2 -q Diff_expression/*_55_39.tsv > Diff_expression/Diff_expression_55_39.tsv", wait=T)
Log_pvalue_39_55 <- read.csv("Diff_expression/Diff_expression_55_39.tsv", 
                             header=T, 
                             sep="\t",
                             quote="",
                             check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_39_55) <- c("Geneid", "Log2chnage_39_55", "padj_39_55")

system("tail -n +2 -q Diff_expression/*_allhydrated_vs_alldry.tsv > Diff_expression/Diff_expression_allhydrated_vs_alldry.tsv", wait=T)
Log_pvalue_allhydrated_vs_alldry <- read.csv("Diff_expression/Diff_expression_allhydrated_vs_alldry.tsv", 
                              header=T, 
                              sep="\t",
                              quote="",
                              check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_allhydrated_vs_alldry) <- c("Geneid", "Log2chnage_allhydrated_vs_alldry", "padj_allhydrated_vs_alldry")

system("tail -n +2 -q Diff_expression/*_earlyhydrated_vs_dry.tsv > Diff_expression/Diff_expression_earlyhydrated_vs_dry.tsv", wait=T)
Log_pvalue_earlyhydrated_vs_dry <- read.csv("Diff_expression/Diff_expression_earlyhydrated_vs_dry.tsv", 
                                             header=T, 
                                             sep="\t",
                                             quote="",
                                             check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_earlyhydrated_vs_dry) <- c("Geneid", "Log2chnage_earlyhydrated_vs_dry", "padj_earlyhydrated_vs_dry")

system("tail -n +2 -q Diff_expression/*_hydrated_vs_earlyhydrated.tsv > Diff_expression/Diff_expression_hydrated_vs_earlyhydrated.tsv", wait=T)
Log_pvalue_hydrated_vs_earlyhydrated <- read.csv("Diff_expression/Diff_expression_hydrated_vs_earlyhydrated.tsv", 
                                            header=T, 
                                            sep="\t",
                                            quote="",
                                            check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_hydrated_vs_earlyhydrated) <- c("Geneid", "Log2chnage_hydrated_vs_earlyhydrated", "padj_hydrated_vs_earlyhydrated")

system("tail -n +2 -q Diff_expression/*_hydrated_vs_dry.tsv > Diff_expression/Diff_expression_hydrated_vs_dry.tsv", wait=T)
Log_pvalue_hydrated_vs_dry <- read.csv("Diff_expression/Diff_expression_hydrated_vs_dry.tsv", 
                                                 header=T, 
                                                 sep="\t",
                                                 quote="",
                                                 check.names = F)[,c(1,3,7)]
colnames(Log_pvalue_hydrated_vs_dry) <- c("Geneid", "Log2chnage_hydrated_vs_dry", "padj_hydrated_vs_dry")

Log_and_pvalues_all <- Reduce(function(x, y) 
  merge(x, y, by="Geneid", all=T), 
  list(Log_pvalue_allhydrated_vs_alldry, Log_pvalue_earlyhydrated_vs_dry, Log_pvalue_hydrated_vs_earlyhydrated, Log_pvalue_hydrated_vs_dry, Log_pvalue_0_0.25, Log_pvalue_0.25_0.5, Log_pvalue_0.5_03, Log_pvalue_03_06, Log_pvalue_06_12, Log_pvalue_12_39, Log_pvalue_39_55)
)

#calculating the smallest p_value. Important to quickly find genes that were significantly up- or down-regulated at lest at some transition.
#na.rm=T is important, because otherwise "NA" will always be the smallest value if it only appears once in a row.

Log_and_pvalues_all <- transform(Log_and_pvalues_all, 
                                 Min_pvalue = pmin(padj_allhydrated_vs_alldry, padj_earlyhydrated_vs_dry, padj_hydrated_vs_earlyhydrated, padj_hydrated_vs_dry, padj_00_0.25, padj_0.25_0.5, padj_0.5_03, padj_03_06, padj_06_12, padj_12_39, padj_39_55, 
                                                   na.rm = T))

write.table(Log_and_pvalues_all, file="Diff_expression/Differential_expression_Log2_and_padj.tsv", 
            sep="\t", 
            quote=F, 
            row.names=F)

All_TPM_annotated <- read.csv("Transcripts_per_bin/Annotated_TPM_perbin_all.tsv", 
                              header=T, 
                              sep="\t",
                              quote="",
                              check.names = F)

ALL_COMBINED <- merge(All_TPM_annotated, Log_and_pvalues_all, 
                      by="Geneid", 
                      all=T)

write.table(ALL_COMBINED, file="Avdat_per_bin_normalized_TPMs_DeSeq2.tsv", 
            sep="\t", 
            quote=F, 
            row.names=F)

#---------------------------------------------------------------------------------------------------------------------------------
# Calculating differentially expressed genes per MAG per timepoint
#---------------------------------------------------------------------------------------------------------------------------------

Avdat_pvalues <- data.table(ALL_COMBINED[,c(1,36,38,40,42,44,46,48,49)])

#replacing missing pvalues with a high number
Avdat_pvalues[is.na(Avdat_pvalues)] <- 100

Diff_expr_perMAG_per_time <- cube(Avdat_pvalues,
                                  lapply(.SD,
                                         function(c)sum(c<=0.05)), 
                                  by="MAG", 
                                  .SDcol=c(2:8))

write.table(Diff_expr_perMAG_per_time, file="Avdat_Diff_expr_perMAG_per_time.tsv", 
            sep="\t", 
            quote=F, 
            row.names=F)