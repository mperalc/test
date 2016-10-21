# Making venn diagrams and mergin DEA tables with peak and across-stages results

#load libraries
library(VennDiagram)



stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7") # 8 stages

currentDate <- Sys.Date() # to save date in name of output files



######### Read in data

#initialize list of data frames to save them inside the loop

DE_stage_specific<-list()
DE_across_stages<-list()

DE_merged<-list()  # for final tables

#### the following might require changes


for(i in stage){
  
  
  DE_stage_specific[[i]] <- read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/conservative_counts/2016-10-17_sig_maxvals_",i,"_diff_expression_results_logFC1.csv",sep=""))
  colnames(DE_stage_specific[[i]]) <- c("ensembl_gene_id", "external_gene_name","gene_biotype", "chromosome_name","logFC","adj.P.Val")
  # get contrast matrix results
  
  
  DE_across_stages[[i]] <- read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/conservative_counts/Contrasts_as_Martijn/2016-10-17_sig_",i,"_diff_expression_maxvals_contrasts_results_logFC1.csv",sep=""))
  
}


  # make list that includes: ensembl_ID, gene, gene_biotype, chr,logFC_peak, logFC_timecourse, adjust.pval_peak, adjust.pval_timecourse, sig_peak/timecourse
  
  try = merge(DE_stage_specific[[i]][c(1:6)],DE_across_stages[[i]][c(1:6)],by=c("ensembl_gene_id","external_gene_name","gene_biotype","chromosome_name"),all=TRUE)
  
  try=try[c(1,2,3,4,5,7,6,8)]
  colnames(try)=c("ensembl_gene_id","gene_name","gene_biotype","chr","logFC_peak","logFC_timecourse","adjust.pval_peak","adjust.pval_timecourse")
  
  # compare 2 pval columns, and make list with presence of sig genes according to NA values. Finally, append list as column 
  
  peak_nas=is.na(try[7])
  timecourse_nas=is.na(try[8])
  
  try[which(!is.na(try[7]) & !is.na(try[8])),9] <- "both"
  try[which(!is.na(try[7]) & is.na(try[8])),9] <- "peak"
  try[which(is.na(try[7]) & !is.na(try[8])),9] <- "timecourse"
  colnames(try)[9]<-"sig_peak/timecourse"
  
  DE_merged[[i]] <- try  # append to list of dataframes
  
  
  
  write.csv(try,quote=F,row.names=F,
            file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/conservative_counts/",
                       currentDate,"_sig_", i,"_diff_expression_maxvals_stage-unique_and_across-stages_logFC1.csv",sep=""))
  
}

# save DE_merged as combined dataframe or as object?





