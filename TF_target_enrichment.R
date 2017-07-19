# This is a test of TF targets enrichment in differential expression analysis for different stages


currentDate <- Sys.Date() # to save date in name of output files

stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")



# test=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/PDX1_both_sig_TIP_results_filtered.txt",
#                     header=T)
#test=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/FOXA2_both_sig_TIP_results_filtered.txt",
 #                   header=T)
# test=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/NKX2_2_both_sig_TIP_results_filtered.txt",
#                      header=T)
# test=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/NKX6_1_both_sig_TIP_results_filtered.txt",
#                       header=T)

# just to save gene lists per stage

PDX1_TIP <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/PDX1_both_sig_TIP_results_filtered.txt",
                                             header=T)
FOXA2_TIP <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/FOXA2_both_sig_TIP_results_filtered.txt",
                        header=T)
NKX6_1_TIP <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/NKX2_2_both_sig_TIP_results_filtered.txt",
                        header=T)
NKX2_2_TIP <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/NKX6_1_both_sig_TIP_results_filtered.txt",
                         header=T)

sig_stages=list()
for(s in stage){
  
  sig_stages[[s]] <- read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages/2017-03-01_sig_",s,"_diff_expression_maxvals_across-stages_results_logFC1.csv",sep=""),header = T) # read in DEA data for each stage
  colnames(sig_stages[[s]])[1]="GeneID"
  sig_stages[[s]]=sig_stages[[s]][ order(sig_stages[[s]]$adj.P.Val,decreasing=F),] # order by adjusted p.value 
  #top 500 genes
  # sig_stages[[s]]=sig_stages[[s]][c(1:500),]
}

# target genes for TF in each stage
target_in_stage_PDX1=list()
target_in_stage_FOXA2=list()
target_in_stage_NKX2_2=list()
target_in_stage_NKX6_1=list()

for(s in stage){
   target_in_stage_PDX1[[s]] = sig_stages[[s]][which(sig_stages[[s]]$GeneID %in% PDX1_TIP$ensembl_gene_id),2]
   # attach TF studied for posterior analyses
   target_in_stage_PDX1[[s]] = c(as.character(target_in_stage_PDX1[[s]]),"PDX1")
   write.table(target_in_stage_PDX1[[s]],
               quote=F,
               row.names=F,col.names = F,
               file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/enrichment_allgenes/gene_lists/PDX1/",
                                                               currentDate,
                                                               "_PDX1_targets_TIP_all_genes_DEA_timecourse_",
                                                                s,
                                                                ".txt",
                                                               sep=""),sep="\t")
   target_in_stage_FOXA2[[s]] = sig_stages[[s]][which(sig_stages[[s]]$GeneID %in% FOXA2_TIP$ensembl_gene_id),2]
   target_in_stage_FOXA2[[s]] = c(as.character(target_in_stage_FOXA2[[s]]),"FOXA2")
   write.table(target_in_stage_FOXA2[[s]],
               quote=F,
               row.names=F,col.names = F,
               file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/enrichment_allgenes/gene_lists/FOXA2/",
                          currentDate,
                          "_FOXA2_targets_TIP_all_genes_DEA_timecourse_",
                          s,
                          ".txt",
                          sep=""),sep="\t")
   target_in_stage_NKX2_2[[s]] = sig_stages[[s]][which(sig_stages[[s]]$GeneID %in% NKX2_2_TIP$ensembl_gene_id),2]
   target_in_stage_NKX2_2[[s]] = c(as.character(target_in_stage_NKX2_2[[s]]),"NKX2_2")
   write.table(target_in_stage_NKX2_2[[s]],
               quote=F,
               row.names=F,col.names = F,
               file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/enrichment_allgenes/gene_lists/NKX2_2/",
                          currentDate,
                          "_NKX2_2_targets_TIP_all_genes_DEA_timecourse_",
                          s,
                          ".txt",
                          sep=""),sep="\t")
   target_in_stage_NKX6_1[[s]] = sig_stages[[s]][which(sig_stages[[s]]$GeneID %in% NKX6_1_TIP$ensembl_gene_id),2]
   target_in_stage_NKX6_1[[s]] = c(as.character(target_in_stage_NKX6_1[[s]]),"NKX6_1")
   
   write.table(target_in_stage_NKX6_1[[s]],
               quote=F,
               row.names=F,col.names = F,
               file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/enrichment_allgenes/gene_lists/NKX6_1/",
                          currentDate,
                          "_NKX6_1_targets_TIP_all_genes_DEA_timecourse_",
                          s,
                          ".txt",
                          sep=""),sep="\t")
}



# hypergeom test to see if there is higher number than expected by chance


library(Homo.sapiens)
library(dplyr)
library(biomaRt)
library(ggplot2)


stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")

load("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_cc.xz")  # 15220 genes and lincRNA


phyperRandom <- function(myGeneList, myGeneSet, genome){
  myRandomGS <- sample( genome,size=length(myGeneSet) )  # take random sample of genes from whole genome, same size as set tested
  myX <- length(which(myGeneList %in% myRandomGS))  # number of T2D in random set
  myM <- length(which(myGeneList %in% genome)) #  T2D genes in my pool
  myN <- length(genome) - length(which(myGeneList %in% genome))  # genes in genome not T2D 
  myK <- length(myGeneSet) # size of set tested
  return(1-phyper(q=myX-1,
                  m=myM, n=myN, k=myK))
}


prob_results=list()
sums_probs=list()
b_test <- list()
overlap <- list()

for(z in stage){
  
  # p-value for enrichment in T2D genes using my differentially expressed genes for each stage
  prob_results[[z]]=1-phyper(q=length(which(test$ensembl_gene_id %in% sig_stages[[z]]$GeneID))-1,
                             m=length(which(test$ensembl_gene_id %in% dge_cc$genes$ensembl_gene_id)),
                             n=length(dge_cc$genes$ensembl_gene_id) - length(which(test$ensembl_gene_id %in% dge_cc$genes$ensembl_gene_id)),
                             k=length(sig_stages[[z]]$GeneID))
  
  # getting amount of T2D genes in each stage and in background pool
  overlap[["total_in_background"]] = length(which(test$ensembl_gene_id %in% dge_cc$genes$ensembl_gene_id))
  
  overlap[[z]] = length(which(test$ensembl_gene_id %in% sig_stages[[z]]$GeneID))
  
  
  
  # same, for random set of genes of same size. random gene set distribution of p-values
  probab_random=numeric()
  
  
  for(i in 1:10000){
    probab_random[i] <- phyperRandom( test$ensembl_gene_id, sig_stages[[z]]$GeneID, dge_cc$genes$ensembl_gene_id )  # get probs from random gene sets, for each diff expr set tested
  }
  
  hist(probab_random)
  
  # empirical p-value 
  
  sums_probs[[z]]$pval=((sum(probab_random<=prob_results[[z]])+1)/(10000+1)) # fraction of probabilities from random draws that are more extreme (smaller) than the one from my tested set
  # +1 in case sum(probab_random<=prob_results[[z]]) ==0
  # that would be the permuted? p-value
  
  # calculate the 95% confidence interval from the binomial distribution
  b_test[[z]]=binom.test(sum(probab_random<=prob_results[[z]]),10000)  # If those fractions are my "successes", then the binomial test can give me a confidence interval
  sums_probs[[z]]$conf.int.low=b_test[[z]]$conf.int[1]
  sums_probs[[z]]$conf.int.up=b_test[[z]]$conf.int[2]
  
  
  print(z)
  
}
df= do.call("rbind", sums_probs)
colnames(df)=c("pval","lowCI","highCI")

overlap_list<- do.call("rbind", overlap)


write.table(df,quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/enrichment_allgenes/",
                                                                      currentDate,
                                                                      "FOXA2_targets_TIP_all_genes_DEA_timecourse_pvals_hyperg_genes_enrichment_10k_permutations.txt",
                                                            sep=""),sep="\t")
write.table(overlap_list,quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/enrichment_allgenes/",
                                                            currentDate,
                                                            "FOXA2_targets_TIP_all_genes_DEA_timecourse_overlap_genes_enrichment_10k_permutations.txt",
                                                            sep=""),sep="\t")
# 
# #### TF target enrichment in correlated genes
# 
# 
# 
# library(Homo.sapiens)
# library(dplyr)
# library(biomaRt)
# library(ggplot2)
# 
# currentDate <- Sys.Date() # to save date in name of output files
# 
# stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")
# 
# 
# load("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_cc.xz")  # 15220 genes and lincRNA
# 
# # gene test data
# 
# PDX1_modules=read.csv2("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/old_counts/conservative_counts/cluster_Agata/WGCNA/WGCNA.connectivity.gene2module_PDX1.csv",
#                      header = F)
# PDX1_modules=as.character(PDX1_modules$V1)
# PDX1_chipseq=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/PDX1_both_sig_TIP_results_filtered.txt",
#                              header=T)
# PDX1_chipseq=as.character(PDX1_chipseq$external_gene_name)
# 
# yes_chip_yes_cor = sum(PDX1_chipseq %in% PDX1_modules)
# 
# # check from here
# yes_chip_not_cor = nrow(dge_cc) - length(setdiff(PDX1_modules,dge_cc$genes$external_gene_name)
# not_chip_yes_cor = sum(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID)
# not_chip_not_cor= nrow(sig_stages[[z]]) - trait_stage
# 
# contingency= matrix(c(trait_stage, trait_background, not_trait_stage, not_trait_background),
#                     nrow = 2,
#                     dimnames =
#                       list(c("stage", "background"),
#                            c("trait", "not trait")))
# 
# test_fisher <- fisher.test(contingency,alternative = "greater")