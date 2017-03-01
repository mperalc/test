#####add info in gene list: which ones are Transcription factors?

currentDate <- Sys.Date() # to save date in name of output files

tf = read.csv('/Users/Marta/Documents/WTCHG/DPhil/Data/TF_list_ensembl.txt',sep=",") # TF list

# add TF yes/no column to files
merged_GWAS_50 <- list()
merged_GWAS_100 <- list()
merged_GWAS_150 <- list()
stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")

for(s in stage){
  
  merged_GWAS_50[[s]] = read.delim(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/GWAS_and_OMIM/50kb/2017-03-01_allinfo_DEA_maxvals_logFC1_50kbs_around_GWAS_loci_",
                                              s,".tsv",sep=""),header=T,check.names=F)
  merged_GWAS_50[[s]][which(merged_GWAS_50[[s]]$ensembl_gene_id %in% tf$Ensembl.Gene.ID),"is_TF"]="yes"
  merged_GWAS_50[[s]][which(!merged_GWAS_50[[s]]$ensembl_gene_id %in% tf$Ensembl.Gene.ID),"is_TF"]="no"
  write.table(merged_GWAS_50[[s]],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/GWAS_and_OMIM/50kb/",
                                                                           currentDate,"_allinfo_DEA_maxvals_logFC1_50kbs_around_GWAS_loci_",s,".tsv",sep=""),sep="\t")
  
  
  merged_GWAS_100[[s]] = read.delim(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/GWAS_and_OMIM/100kb/2017-03-01_allinfo_DEA_maxvals_logFC1_100kbs_around_GWAS_loci_",
                                               s,".tsv",sep=""),header=T,check.names=F)
  merged_GWAS_100[[s]][which(merged_GWAS_100[[s]]$ensembl_gene_id %in% tf$Ensembl.Gene.ID),"is_TF"]="yes"
  merged_GWAS_100[[s]][which(!merged_GWAS_100[[s]]$ensembl_gene_id %in% tf$Ensembl.Gene.ID),"is_TF"]="no"
  write.table(merged_GWAS_100[[s]],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/GWAS_and_OMIM/100kb/",
                                                                            currentDate,"_allinfo_DEA_maxvals_logFC1_100kbs_around_GWAS_loci_",s,".tsv",sep=""),sep="\t")
  
  
  merged_GWAS_150[[s]] = read.delim(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/GWAS_and_OMIM/150kb/2017-03-01_allinfo_DEA_maxvals_logFC1_150kbs_around_GWAS_loci_",
                                               s,".tsv",sep=""),header=T,check.names=F)
  merged_GWAS_150[[s]][which(merged_GWAS_150[[s]]$ensembl_gene_id %in% tf$Ensembl.Gene.ID),"is_TF"]="yes"
  merged_GWAS_150[[s]][which(!merged_GWAS_150[[s]]$ensembl_gene_id %in% tf$Ensembl.Gene.ID),"is_TF"]="no"
  write.table(merged_GWAS_150[[s]],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/GWAS_and_OMIM/150kb/",
                                                                            currentDate,"_allinfo_DEA_maxvals_logFC1_150kbs_around_GWAS_loci_",s,".tsv",sep=""),sep="\t")
}






