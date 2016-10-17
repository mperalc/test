# Making venn diagrams

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
  

# Save venn plot for every stage
  
  
  for(i in stage){
    
  
# need to get my results into three variables:
  # stage-unique= x 
  #combined = y
  # across-stages= z
#by distinguishing which ones are shared, and which ones are found in DE_stage_specific and not in the other database and viceversa.
    
    x=setdiff(DE_stage_specific[[i]][[1]],DE_across_stages[[i]][[1]])
    y=intersect(DE_stage_specific[[i]][[1]],DE_across_stages[[i]][[1]])
    z=setdiff(DE_across_stages[[i]][[1]],DE_stage_specific[[i]][[1]])
    
    
    #merge into list
    list=c(x,y,z)
    # stage_specific goes from 1 to lenght(x)+length(y)
    # across_stages goes from lenght(x)+1 to length(list)
    setwd("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/venn_plots_contrasts_normalDE")  
    venn.plot <- venn.diagram(x=list(A=1:(length(x)+length(y)),B=(length(x)+1):length(list)),
                              filename = paste("Venn_",i,"_stage_DE_stage-specific_and_across-stages_logFC1_",currentDate,".tiff",sep=""),
                              lwd=c(0.8),
                              fill=c("darkmagenta", "darkblue"),
                              alpha=c(0.5,0.5), 
                              cex = 1, 
                              cat.fontface=2, 
                              cat.dist = c(0.06, 0.06),
                              cat.pos = c(-90, 90),
                              cat.cex=1,
                              # cat.default.pos=c("outer"),
                              category.names=c("stage specific", "across stages"),
                              # col = "black",
                              # 
                              scaled = TRUE,
                              # ext.text = TRUE,
                              # ext.line.lwd = 1,
                              ext.dist = 0.01,
                              ext.length = 0.9,
                              #  #ext.pos = -4,
                              # inverted = TRUE,
                              margin=0.1,
                              
                              main = paste(i,"stage",sep=" "),
                              sub = "Significantly DE genes assigned to this stage for each method",
                              main.fontface=2,
                              main.cex = 2,
                              main.pos = c(0.5,1.06),
                              sub.cex = 1,
                              sub.pos=  c(0.5,1)
    )
    
  

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
  
  
  
  
  
