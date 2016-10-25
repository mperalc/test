# annotating genes around (150kb) GWAS loci
# diabetes and glycemic traits

library("biomaRt")

currentDate <- Sys.Date() # to save date in name of output files


############ read in GWAS loci #######################

T2D <- read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/T2D_GWAS_loci_061115", header = F, sep = "\t")
colnames(T2D) = c("SNP","GWAS_loci","location")
glycemic <- read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Glycemic_GWAS_loci_061115", header = F, sep = "\t")
colnames(glycemic) = c("SNP","GWAS_loci","location")

shared_SNPs = T2D[match(glycemic$SNP,T2D$SNP),1] # list of SNPs in T2D list that are shared by glycemic list
shared_SNPs=shared_SNPs[which(!is.na(shared_SNPs))]  #eliminate NAs
# to be marked with "both"

#eliminate duplicated SNPs in glycemic df
glycemic = glycemic[!glycemic$SNP %in% shared_SNPs,]

GWAS = rbind(T2D,glycemic) # merging
trait=c(rep("T2D",nrow(T2D)),rep("glycemic",nrow(glycemic))) # labeling GWAS trait (T2D or glycemic)
GWAS$trait=trait

GWAS[which(GWAS$SNP %in% shared_SNPs),4] = "both"  # marking SNPs shared by both traits

GWAS = GWAS[which(!is.na(GWAS$location)),]# remove SNPs with NA in location (sex chromosomes)

GWAS$location=gsub("\\chr*","",GWAS$location) #remove "chr" in every element in location



#########add and substract n kb to get chr:xxx-n:xxx+n format #################

kbs <- list()

distance = c("50","100","150") # list of distances /1000 (kb)

for(d in distance){
  
  kbs[[d]] <- GWAS
  minus_d =  as.numeric(gsub(".*:","",kbs[[d]]$location))-as.numeric(d)*1000          #vector of -distance
  plus_d =  as.numeric(gsub(".*:","",kbs[[d]]$location))+as.numeric(d)*1000   #vector of +distance
  
  kbs[[d]]$location = paste(gsub("\\:.*","",kbs[[d]]$location),minus_d,plus_d,sep=":") #paste everything
  rm(minus_d,plus_d)
  write.table(kbs[[d]],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/",
                                                       currentDate,"_location_",d,"kbs_around_GWAS_loci_glycemic_and_T2D.tsv",sep=""),sep="\t")
}



### use those kbs lists to do query

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
# the lists of GWAS hits I was given use genome coordinates from GRCh37 (not GRCh38, which is the latest release). Keep that in mind.

listMarts(ensembl)

filters = listFilters(ensembl)


results <- list()

for(d in distance){
  filterlist <- as.list(kbs[[d]]$location)
  results[[d]] <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name', "start_position", "end_position"),
              filters = c("chromosomal_region","biotype"),values = list(chromosomal_region=filterlist, biotype = c("protein_coding","lincRNA")), mart = ensembl)
  
  write.table(results[[d]],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/",
                                                      currentDate,"_genes_and_lincRNA",d,"kbs_around_GWAS_loci_glycemic_and_T2D.tsv",sep=""),sep="\t")
}


## look for T2D associated gene enrichment with Fisher


# import dge object

load(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_cc.xz",verbose=TRUE)  #loading the dge object for conservative counts


# load files with DE results
sig_stages=list()
stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")

for(s in stage){
  # sig_stages[[i]]= read.csv(paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/conservative_counts/2016-10-17_sig_maxvals_",i,"_diff_expression_results_logFC1.csv",sep=""))
  # colnames(sig_stages[[i]]) <- c("ensembl_gene_id", "external_gene_name","gene_biotype", "chromosome_name","logFC","adj.P.Val")
  # 
  sig_stages[[s]]= read.delim(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/conservative_counts/comparison_with_and_without_contrasts/OMIM/2016-10-24_sig_", 
                                       s,"_DEA_maxvals_stage-unique_and_across-stages_logFC1_OMIM.tsv",sep=""),header=T,check.names=F)
  sig_stages[[s]]$ensembl_gene_id = as.character(sig_stages[[s]]$ensembl_gene_id) # for posterior comparison, ensembl ids here and in result table need to be characters
}


s="iPSC"
#iPSC

iPSC_t2d_50 <- fisher.test(data.frame(GWAS_loci=c(length(which(unique(sig_stages[[s]]$ensembl_gene_id) %in% results$`50`$ensembl_gene_id)),length(which(rownames(dge_cc) %in%  results$`50`$ensembl_gene_id))),
                                      not_GWAS_loci=c(length(unique(sig_stages[[s]]$ensembl_gene_id)) - length(which(unique(sig_stages[[s]]$ensembl_gene_id) %in% results$`50`$ensembl_gene_id)),dim(dge_cc)[1] - length(which(rownames(dge_cc) %in%  results$`50`$ensembl_gene_id)))))

iPSC_t2d_100 <- fisher.test(data.frame(GWAS_loci=c(length(which(unique(sig_stages[[s]]$ensembl_gene_id) %in% results$`100`$ensembl_gene_id)),length(which(rownames(dge_cc) %in%  results$`100`$ensembl_gene_id))),
                                       not_GWAS_loci=c(length(unique(sig_stages[[s]]$ensembl_gene_id)) - length(which(unique(sig_stages[[s]]$ensembl_gene_id) %in% results$`100`$ensembl_gene_id)),dim(dge_cc)[1] - length(which(rownames(dge_cc) %in%  results$`100`$ensembl_gene_id)))))


iPSC_t2d_150 <- fisher.test(data.frame(GWAS_loci=c(length(which(unique(sig_stages[[s]]$ensembl_gene_id) %in% results$`150`$ensembl_gene_id)),length(which(rownames(dge_cc) %in%  results$`150`$ensembl_gene_id))),
                                       not_GWAS_loci=c(length(unique(sig_stages[[s]]$ensembl_gene_id)) - length(which(unique(sig_stages[[s]]$ensembl_gene_id) %in% results$`150`$ensembl_gene_id)),dim(dge_cc)[1] - length(which(rownames(dge_cc) %in%  results$`150`$ensembl_gene_id)))))

s="EN7"

EN7_t2d_50 <- fisher.test(data.frame(GWAS_loci=c(length(which(unique(sig_stages[[s]]$ensembl_gene_id) %in% results$`50`$ensembl_gene_id)),length(which(rownames(dge_cc) %in%  results$`50`$ensembl_gene_id))),
                                      not_GWAS_loci=c(length(unique(sig_stages[[s]]$ensembl_gene_id)) - length(which(unique(sig_stages[[s]]$ensembl_gene_id) %in% results$`50`$ensembl_gene_id)),dim(dge_cc)[1] - length(which(rownames(dge_cc) %in%  results$`50`$ensembl_gene_id)))))

EN7_t2d_100 <- fisher.test(data.frame(GWAS_loci=c(length(which(unique(sig_stages[[s]]$ensembl_gene_id) %in% results$`100`$ensembl_gene_id)),length(which(rownames(dge_cc) %in%  results$`100`$ensembl_gene_id))),
                                       not_GWAS_loci=c(length(unique(sig_stages[[s]]$ensembl_gene_id)) - length(which(unique(sig_stages[[s]]$ensembl_gene_id) %in% results$`100`$ensembl_gene_id)),dim(dge_cc)[1] - length(which(rownames(dge_cc) %in%  results$`100`$ensembl_gene_id)))))


EN7_t2d_150 <- fisher.test(data.frame(GWAS_loci=c(length(which(unique(sig_stages[[s]]$ensembl_gene_id) %in% results$`150`$ensembl_gene_id)),length(which(rownames(dge_cc) %in%  results$`150`$ensembl_gene_id))),
                                       not_GWAS_loci=c(length(unique(sig_stages[[s]]$ensembl_gene_id)) - length(which(unique(sig_stages[[s]]$ensembl_gene_id) %in% results$`150`$ensembl_gene_id)),dim(dge_cc)[1] - length(which(rownames(dge_cc) %in%  results$`150`$ensembl_gene_id)))))


# do this for peak/timecourse DEA separately