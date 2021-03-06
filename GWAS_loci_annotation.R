# annotating genes around (150kb) GWAS loci
# diabetes and glycemic traits

library("biomaRt")
library("IRanges")

currentDate <- Sys.Date() # to save date in name of output files


distance = c("50","100","150") # list of distances /1000 (kb)

############ read in GWAS loci #######################
# 
# T2D <- read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/T2D_GWAS_loci_061115", header = F, sep = "\t")
# colnames(T2D) = c("SNP","GWAS_loci","location")
# glycemic <- read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Glycemic_GWAS_loci_061115", header = F, sep = "\t")
# colnames(glycemic) = c("SNP","GWAS_loci","location")
# 
# shared_SNPs = T2D[match(glycemic$SNP,T2D$SNP),1] # list of SNPs in T2D list that are shared by glycemic list
# shared_SNPs=shared_SNPs[which(!is.na(shared_SNPs))]  #eliminate NAs
# # to be marked with "both"
# 
# #eliminate duplicated SNPs in glycemic df
# glycemic = glycemic[!glycemic$SNP %in% shared_SNPs,]
# 
# GWAS = rbind(T2D,glycemic) # merging
# trait=c(rep("T2D",nrow(T2D)),rep("glycemic",nrow(glycemic))) # labeling GWAS trait (T2D or glycemic)
# GWAS$trait=trait
# 
# GWAS[which(GWAS$SNP %in% shared_SNPs),4] = "both"  # marking SNPs shared by both traits
# 
# GWAS = GWAS[which(!is.na(GWAS$location)),]# remove SNPs with NA in location (sex chromosomes)
# 
# GWAS$location=gsub("\\chr*","",GWAS$location) #remove "chr" in every element in location
# 
# 
# 
# #########add and substract n kb to get chr:xxx-n:xxx+n format #################
# 
# kbs <- list()

# for(d in distance){
#   
#   kbs[[d]] <- GWAS
#   minus_d =  as.numeric(gsub(".*:","",kbs[[d]]$location))-as.numeric(d)*1000          #vector of -distance
#   plus_d =  as.numeric(gsub(".*:","",kbs[[d]]$location))+as.numeric(d)*1000   #vector of +distance
#   
#   kbs[[d]]$location = paste(gsub("\\:.*","",kbs[[d]]$location),minus_d,plus_d,sep=":") #paste everything
#   rm(minus_d,plus_d)
#   write.table(kbs[[d]],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/",
#                                                        currentDate,"_location_",d,"kbs_around_GWAS_loci_glycemic_and_T2D.tsv",sep=""),sep="\t")
# }
# 
# 
# 
# ### use those kbs lists to do query
# 
# ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
# # the lists of GWAS hits I was given use genome coordinates from GRCh37 (not GRCh38, which is the latest release). Keep that in mind.
# 
# listMarts(ensembl)
# 
# filters = listFilters(ensembl)
# 
# 
# GWAS_list <- list()
# 
# 
# for(d in distance){
#   filterlist <- as.list(kbs[[d]]$location)
#   GWAS_list[[d]] <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name', "start_position", "end_position"),
#               filters = c("chromosomal_region","biotype"),
#               values = list(chromosomal_region=filterlist, biotype = c("protein_coding","lincRNA")),
#               mart = ensembl)
# 
#   GWAS_list[[d]] <- GWAS_list[[d]][order(GWAS_list[[d]]$chromosome_name),]
#   
# }
# 
# 
# 
# 
# 
# ############ merge genes and GWAS loci they come from
# kbs_simple <- list()
# merged <- list()
# 
# 
# for(d in distance){
#   kbs_simple[[d]] <- kbs[[d]][,c(1,2,4)]
#   kbs_simple[[d]]$chromosome_name= as.numeric(gsub("(.*?):.*", "\\1",kbs[[d]]$location))
#   minus_d= gsub("^[^:]*:","",kbs[[d]]$location)             # remove everything up to 1st colon
#   kbs_simple[[d]]$start_position =  as.numeric(gsub("(.*?):.*", "\\1",minus_d))    # remove everything after colon
#   rm(minus_d)
#   kbs_simple[[d]]$end_position =  as.numeric(gsub(".*:","",kbs[[d]]$location))    # remove everything before last colon
# 
#   # creating IRanges objects for easy comparison 
#   kbs_ranges=split(IRanges(kbs_simple[[d]]$start_position,
#                            kbs_simple[[d]]$end_position),
#                    kbs_simple[[d]]$chromosome_name) # range for kbs
#   
#   GWAS_ranges=split(IRanges(GWAS_list[[d]]$start_position, 
#                             GWAS_list[[d]]$end_position,
#                             names = GWAS_list[[d]]$external_gene_name),
#                     GWAS_list[[d]]$chromosome_name) # range for GWAS
#   
#   for(c in 1:20){
#     
#     overlap=findOverlaps(GWAS_ranges[[c]],kbs_ranges[[c]]) # for every chromosome, relates indices in GWAS_ranges with indices in kbs_ranges lists
#  
#     
#     check_query = GWAS_ranges[[c]][overlap@from]@NAMES # genes in location names from overlap indices
#     check_subject = kbs_simple[[d]][match(kbs_ranges[[c]][overlap@to]@start,kbs_simple[[d]]$start_position),]#gets positions in kbs table from start position in subject (got from overlap matches) 
#     check_subject$genes_in_loci=check_query #merges genes in loci (query table) with GWAS_loci table
#     
#     merged[[d]] <- rbind(merged[[d]],check_subject)
#     
#     
#   
#   }
#   merged[[d]]$ensembl_gene_id <- GWAS_list[[d]][match(merged[[d]]$genes_in_loci,GWAS_list[[d]]$external_gene_name),1] ###append gene id of every gene
#   
#   write.table(merged[[d]],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/",
#                                                             currentDate,"_genes_and_lincRNA",d,"kbs_around_GWAS_loci_glycemic_and_T2D.tsv",sep=""),sep="\t")
#   
# }

# read in table with GWAS gene info:
merged <- list()

for(d in distance){
  merged[[d]] <- read.table(file = paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/2016-11-02_genes_and_lincRNA",d,"kbs_around_GWAS_loci_glycemic_and_T2D.tsv",sep=""), header = T, sep = "\t")

}
# write GWAS info in DEA results table

# load files with DEA results
sig_stages=list()
stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")

for(s in stage){
  
  sig_stages[[s]]= read.delim(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/OMIM/2017-03-01_sig_", 
                                       s,"_DEA_maxvals_stage-unique_and_across-stages_logFC1_OMIM.tsv",sep=""),header=T,check.names=F)
  sig_stages[[s]]$ensembl_gene_id = as.character(sig_stages[[s]]$ensembl_gene_id) # for posterior comparison, ensembl ids here and in result table need to be characters
}

# merge and save info in file
merged_GWAS_50 <- list()
merged_GWAS_100 <- list()
merged_GWAS_150 <- list()

  for(s in stage){
    merged_GWAS_50[[s]] = merge(sig_stages[[s]],merged$`50`, by="ensembl_gene_id",all.x=TRUE) 
    merged_GWAS_100[[s]] = merge(sig_stages[[s]],merged$`100`, by="ensembl_gene_id",all.x=TRUE) 
    merged_GWAS_150[[s]] = merge(sig_stages[[s]],merged$`150`, by="ensembl_gene_id",all.x=TRUE) 
  }


#write down file with OMIM and GWAS info
  
for(s in stage){
  
  write.table(merged_GWAS_50[[s]][,c(1:15)],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/GWAS_and_OMIM/50kb/",
                                                         currentDate,"_allinfo_DEA_maxvals_logFC1_50kbs_around_GWAS_loci_",s,".tsv",sep=""),sep="\t")
  
  write.table(merged_GWAS_100[[s]][,c(1:15)],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/GWAS_and_OMIM/100kb/",
                                                                           currentDate,"_allinfo_DEA_maxvals_logFC1_100kbs_around_GWAS_loci_",s,".tsv",sep=""),sep="\t")
  write.table(merged_GWAS_150[[s]][,c(1:15)],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/GWAS_and_OMIM/150kb/",
                                                                           currentDate,"_allinfo_DEA_maxvals_logFC1_150kbs_around_GWAS_loci_",s,".tsv",sep=""),sep="\t")
}
# 
# table_150=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/2016-11-02_genes_and_lincRNA150kbs_around_GWAS_loci_glycemic_and_T2D.tsv",sep=""),sep="\t",header=T)
# table_50=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/2016-11-02_genes_and_lincRNA50kbs_around_GWAS_loci_glycemic_and_T2D.tsv",sep=""),sep="\t",header=T)
# table_100=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/2016-11-02_genes_and_lincRNA100kbs_around_GWAS_loci_glycemic_and_T2D.tsv",sep=""),sep="\t",header=T)
# 
# 
# 
# ## look for T2D associated gene enrichment with Fisher
# 
# stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")
# 
# # import dge object
# 
# load(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_cc.xz",verbose=TRUE)  #loading the dge object for conservative counts
# 
# 
# # split peak and timecourse results for test
# peak <-  lapply(sig_stages, subset, sig_peak.timecourse=="peak" | sig_peak.timecourse=="both" )
# timecourse <- lapply(sig_stages, subset, sig_peak.timecourse=="timecourse" | sig_peak.timecourse=="both" )
# 
# 
# 
# 
# # GWAS peak & timecourse: all stages
# GWAS_50_peak <- list()
# GWAS_50_timecourse <- list()
# GWAS_100_peak <- list()
# GWAS_100_timecourse <- list()
# GWAS_150_peak <- list()
# GWAS_150_timecourse <- list()
# 
# for(s in stage){
#  
#   GWAS_50_peak[[s]] <- fisher.test(data.frame(GWAS_loci=c(length(which(unique(peak[[s]]$ensembl_gene_id) %in% GWAS_list$`50`$ensembl_gene_id)),length(which(rownames(dge_cc) %in%  GWAS_list$`50`$ensembl_gene_id))),
#                                       not_GWAS_loci=c(length(unique(peak[[s]]$ensembl_gene_id)) - length(which(unique(peak[[s]]$ensembl_gene_id) %in% GWAS_list$`50`$ensembl_gene_id)),dim(dge_cc)[1] - length(which(rownames(dge_cc) %in%  GWAS_list$`50`$ensembl_gene_id)))))
#   
#   GWAS_50_timecourse[[s]] <- fisher.test(data.frame(GWAS_loci=c(length(which(unique(timecourse[[s]]$ensembl_gene_id) %in% GWAS_list$`50`$ensembl_gene_id)),length(which(rownames(dge_cc) %in%  GWAS_list$`50`$ensembl_gene_id))),
#                                          not_GWAS_loci=c(length(unique(timecourse[[s]]$ensembl_gene_id)) - length(which(unique(timecourse[[s]]$ensembl_gene_id) %in% GWAS_list$`50`$ensembl_gene_id)),dim(dge_cc)[1] - length(which(rownames(dge_cc) %in%  GWAS_list$`50`$ensembl_gene_id)))))
#   
#   
#   GWAS_100_peak[[s]] <- fisher.test(data.frame(GWAS_loci=c(length(which(unique(peak[[s]]$ensembl_gene_id) %in% GWAS_list$`100`$ensembl_gene_id)),length(which(rownames(dge_cc) %in%  GWAS_list$`100`$ensembl_gene_id))),
#                                               not_GWAS_loci=c(length(unique(peak[[s]]$ensembl_gene_id)) - length(which(unique(peak[[s]]$ensembl_gene_id) %in% GWAS_list$`100`$ensembl_gene_id)),dim(dge_cc)[1] - length(which(rownames(dge_cc) %in%  GWAS_list$`100`$ensembl_gene_id)))))
#   
#   GWAS_100_timecourse[[s]] <- fisher.test(data.frame(GWAS_loci=c(length(which(unique(timecourse[[s]]$ensembl_gene_id) %in% GWAS_list$`100`$ensembl_gene_id)),length(which(rownames(dge_cc) %in%  GWAS_list$`100`$ensembl_gene_id))),
#                                                     not_GWAS_loci=c(length(unique(timecourse[[s]]$ensembl_gene_id)) - length(which(unique(timecourse[[s]]$ensembl_gene_id) %in% GWAS_list$`100`$ensembl_gene_id)),dim(dge_cc)[1] - length(which(rownames(dge_cc) %in%  GWAS_list$`100`$ensembl_gene_id)))))
#   
#   
#   GWAS_150_peak[[s]] <- fisher.test(data.frame(GWAS_loci=c(length(which(unique(peak[[s]]$ensembl_gene_id) %in% GWAS_list$`150`$ensembl_gene_id)),length(which(rownames(dge_cc) %in%  GWAS_list$`150`$ensembl_gene_id))),
#                                               not_GWAS_loci=c(length(unique(peak[[s]]$ensembl_gene_id)) - length(which(unique(peak[[s]]$ensembl_gene_id) %in% GWAS_list$`150`$ensembl_gene_id)),dim(dge_cc)[1] - length(which(rownames(dge_cc) %in%  GWAS_list$`150`$ensembl_gene_id)))))
#   
#   GWAS_150_timecourse[[s]] <- fisher.test(data.frame(GWAS_loci=c(length(which(unique(timecourse[[s]]$ensembl_gene_id) %in% GWAS_list$`150`$ensembl_gene_id)),length(which(rownames(dge_cc) %in%  GWAS_list$`150`$ensembl_gene_id))),
#                                                     not_GWAS_loci=c(length(unique(timecourse[[s]]$ensembl_gene_id)) - length(which(unique(timecourse[[s]]$ensembl_gene_id) %in% GWAS_list$`150`$ensembl_gene_id)),dim(dge_cc)[1] - length(which(rownames(dge_cc) %in%  GWAS_list$`150`$ensembl_gene_id)))))
#   
# }
# 
# # save 2x2 tables and p-values
