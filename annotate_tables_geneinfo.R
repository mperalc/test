# This is part of the paper pipeline

# annotate gene tables
# with whatever gene info we want, specially T2D and T1D genes
# also TFs
library(biomaRt)

currentDate <- Sys.Date() # to save date in name of output files

stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")

sig_stages=list()

# TF targets from Pasquali 2013 hg18


PDX1=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/PDX1_both_sig_TIP_results_filtered.txt",
                header=T)
FOXA2=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/FOXA2_both_sig_TIP_results_filtered.txt",
                header=T)
NKX2_2=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/NKX2_2_both_sig_TIP_results_filtered.txt",
                header=T)
NKX6_1=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/Pasquali_2013/original_hg18/NKX6_1_both_sig_TIP_results_filtered.txt",
                header=T)

#  monogenic diabetes genes Fuchsberger 2016:
monogenic <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Monogenic_diabetes/Monogenic_primary_Fuchsberger_16.txt")
colnames(monogenic)="gene"

#  T2D genes 0kb credible regions DIAGRAM 2017:
T2D <-  read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_T2D/2017-02-22annotated_genes_in_credible_regions_plusminus_0_kb.txt",header=T)

# FG genes 0kb credible regions ENGAGE 2017:

FG <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_FG/2017-05-03FG_annotated_genes_in_credible_regions_plusminus_0_kb.txt",header=T)

# transcription factors

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 


filterlist=c("GO:0003700")
is_tf_GO <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype',
                                     'chromosome_name', "start_position", "end_position","go_id","name_1006"),  
                      filters = c("go"),values=filterlist, mart = ensembl)
is_tf_GO=is_tf_GO[which(is_tf_GO$go_id=="GO:0003700"),]

for(s in stage){
  # read in DEA data for each stage
  sig_stages[[s]] <- read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages/2017-07-03_sig_",
                                         s,"_diff_expression_maxvals_across-stages_results_logFC1_annotated.csv",sep=""),header = T) 

  sig_stages[[s]]=sig_stages[[s]][ order(sig_stages[[s]]$adj.P.Val,decreasing=F),] # order by adjusted p.value of "timecourse" (across-stages) method
  
  sig_stages[[s]]$T2D_credset <-  ifelse(  sig_stages[[s]]$external_gene_name %in% T2D$external_gene_name, 'yes' , 'no' ) 
  sig_stages[[s]]$monogenic_primary <-  ifelse(sig_stages[[s]]$external_gene_name %in% monogenic$gene, 'yes' , 'no' ) 
  sig_stages[[s]]$FG_credset<-  ifelse(sig_stages[[s]]$external_gene_name %in% FG$external_gene_name, 'yes' , 'no' ) 
  sig_stages[[s]]$TF_GO<-  ifelse(sig_stages[[s]]$external_gene_name %in% is_tf_GO$external_gene_name, 'yes' , 'no' ) 
  sig_stages[[s]]$PDX1_target_Pasquali <- ifelse(sig_stages[[s]]$external_gene_name %in% PDX1$external_gene_name, 'yes' , 'no' ) 
  sig_stages[[s]]$FOXA2_target_Pasquali <- ifelse(sig_stages[[s]]$external_gene_name %in% FOXA2$external_gene_name, 'yes' , 'no' ) 
  sig_stages[[s]]$NKX2_2_target_Pasquali <- ifelse(sig_stages[[s]]$external_gene_name %in% NKX2_2$external_gene_name, 'yes' , 'no' ) 
  sig_stages[[s]]$NKX6_1_target_Pasquali <- ifelse(sig_stages[[s]]$external_gene_name %in% NKX6_1$external_gene_name, 'yes' , 'no' ) 
  
  
   write.csv(sig_stages[[s]],file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages/2017-07-03_sig_",s,
                                       "_diff_expression_maxvals_across-stages_results_logFC1_annotated.csv",sep=""),
            row.names = F,quote = F)
  }
