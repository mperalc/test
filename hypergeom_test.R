## Hypergeometric test for datasets: are they enriched in T2D/FG genes?

# Takes in:
  # A list of SNPs in LD with credible sets for DIAGRAM (T2D)
      # I then need to take the SNPs coordinates of the extremes of each regions, and take out all genes inside (ENSEMBL)
  # Datasets to test for enrichment in the above
  # The background dataset to create random sets of genes. It would be the list of genes from the initial RNA-seq list (before filtering)

# Then, I test how my datasets are enriched in T2D genes compared to random sets of genes of the same size. 

library(Homo.sapiens)
library(dplyr)
library(biomaRt)
library(ggplot2)

currentDate <- Sys.Date() # to save date in name of output files

stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")

trait= "FG"    # T2D or fasting glucose (FG)?
### files

#### background "genome" data


# Includes all sequenced genes, including those that have 0 counts:
# background_genes=as.data.frame(read_tsv("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/other_counts/27.07.2016.Differentiation_v2.gene.counts.tsv"))
# background_genes$GeneID <- gsub("(ENSG[0-9]+).*", "//1", background_genes$GeneID) # cut everything after the dot
# background_genes=background_genes$GeneID   # 63568 genes and lincRNA
# 
# # take out those not in ensembl, those in the sex chromosomes and 
# ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
# # I use genome coordinates from GRCh37 (not GRCh38, which is the latest release). Keep that in mind.
# # that's why I don't want to have ids that have dissappeared in the latest build
# 
# all_ensembl_info_cc <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name'),
#                              filters = 'ensembl_gene_id', values = background_genes, mart = ensembl)
# 
# background_genes <- background_genes[all_ensembl_info_cc$ensembl_gene_id] #keep genes in counts whose geneID still exists in Ensembl

# just data used for the Differential expression analysis

load("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_cc.xz")  # 15220 genes and lincRNA

# gene test data

distance = c("50","100","200","500") # list of distances in kb

test <- list()
test[["0"]]=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                  trait,
                                  # "/2017-02-22_annotated_genes_in_credible_regions.txt",sep=""),header=T)  # T2D
                                  "/FG_genes_credible_set_0kb.txt",sep=""),header=T)
colnames(test[["0"]])[1]=c("GeneID")
test[["0"]]$GeneID=as.character(test[["0"]]$GeneID)                    

for(d in distance){
  test[[d]]=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                  trait,
                                  "/2017-02-22annotated_genes_in_credible_regions_plusminus_",
                                  d,"_kb.txt",sep=""),header=T)
  colnames(test[[d]])[1]=c("GeneID")
  test[[d]]$GeneID=as.character(test[[d]]$GeneID)
}


# files with DE results (across-stages)

# # all genes, get list of ensembl gene ids

sig_stages=list()
for(s in stage){
 
  sig_stages[[s]] <- read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages/2017-03-01_sig_",s,"_diff_expression_maxvals_across-stages_results_logFC1.csv",sep=""),header = T) # read in DEA data for each stage
  colnames(sig_stages[[s]])[1]="GeneID"
  sig_stages[[s]]=sig_stages[[s]][ order(sig_stages[[s]]$adj.P.Val,decreasing=F),] # order by adjusted p.value 
  # top 500 genes
 # sig_stages[[s]]=sig_stages[[s]][c(1:500),]
  }

##################
# hypergeometric test with permutations

# Total balls in the urn: background dataset (all RNA-seq genes)
# x --> number of white balls drawn without replacement from an urn which contains both black and white balls (T2D genes in my tested dataset)
# m --> White balls in the urn: characteristic I'm testing for enrichment (T2D genes)
# n --> Black balls in the urn: total-white balls (non-T2D genes)
# k --> Number of balls drawn from the urn: depends on the size of the dataset tested 

# I want the p-value of getting (number of T2D genes in my DE dataset) or more white balls in a sample 
# of size (number of genes in my DE dataset) from an urn with (number of T2D genes) white balls and 
# (total genes RNAseq-T2D genes in that set) black balls.

# 1-phyper(q, m, n, k)

# q in this case is x-1 (diff exp dataset -1), the probability of x or more. It is a quantile defined as the smallest value x such that F(x) ≥ p, 
# where F is the distribution function.
# 
# 
# dhyperRandom <- function(myGeneList, myGeneSet, genome){
#   myRandomGS <- sample( genome,size=length(myGeneSet) )  # take random sample of genes from whole genome, same size as set tested
#   myX <- length(which(myGeneList %in% myRandomGS))  # number of T2D in random set
#   myM <- length(which(myGeneList %in% genome)) #  T2D genes in my pool
#   myN <- length(genome) - length(which(myGeneList %in% genome))  # genes in genome not T2D 
#   myK <- length(myGeneSet) # size of set tested
#   return(1-phyper(q=myX-1, m=myM, n=myN, k=myK))
#   #return(dhyper(x=myX, m=myM, n=myN, k=myK))
# }
# 
# random_dist=numeric()
# randomGS_dist <- function(myGeneList,myGeneSet,genome){
#   myRandomGS= sample(genome,size=length(myGeneSet))
#   return(length(which(myGeneList %in% myRandomGS)) ) 
# }
# 
# for(i in 1:1000){
#   random_dist[i] <- randomGS_dist( T2D$GeneID, sig_stages$EN7$GeneID, dge_cc$genes$ensembl_gene_id )  # get p-values from random gene sets, for each diff expr set tested
# }
# 
# hist(random_dist)
# results_EN7_dist= length(which(T2D$GeneID %in% sig_stages$EN7$GeneID))
# 
# sum(random_dist>=results_EN7_dist)/1000

# p vals directly 
# 
# 
# pvalue_random=numeric()
# 
# for(i in 1:1000){
#   pvalue_random[i] <- dhyperRandom( T2D$GeneID, sig_stages$EN7$GeneID, dge_cc$genes$ensembl_gene_id )  # get p-values from random gene sets, for each diff expr set tested
# }
# 
# hist(pvalue_random)
# # now I have a distribution of p-values, and I want to test whether that p-value is in the 5% most extreme tail
# 
# # observed p-values:
# pval_results_EN7=1-phyper(q=length(which(T2D$GeneID %in% sig_stages$EN7$GeneID))-1,
#                      m=length(which(T2D$GeneID %in% dge_cc$genes$ensembl_gene_id)),
#                      n=length(dge_cc$genes$ensembl_gene_id) - length(which(T2D$GeneID %in% dge_cc$genes$ensembl_gene_id)),
#                      k=length(sig_stages$EN7$GeneID))
# # 
# # results_EN7=dhyper(x=length(which(T2D$GeneID %in% sig_stages$EN7$GeneID)),
# #                                           m=length(which(T2D$GeneID %in% dge_cc$genes$ensembl_gene_id)),
# #                                           n=length(dge_cc$genes$ensembl_gene_id) - length(which(T2D$GeneID %in% dge_cc$genes$ensembl_gene_id)),
# #                                           k=length(sig_stages$EN7$GeneID))
# 
# sum(pvalue_random<=pval_results_EN7)/1000



# repeat for every diff expr gene set result

# comparison of probabilities

df_list <- list()
distance = c("0","50","100","200","500") # list of distances in kb


for(d in distance){
  # dhyperRandom <- function(myGeneList, myGeneSet, genome){
  #   myRandomGS <- sample( genome,size=length(myGeneSet) )  # take random sample of genes from whole genome, same size as set tested
  #   myX <- length(which(myGeneList %in% myRandomGS))  # number of T2D in random set
  #   myM <- length(which(myGeneList %in% genome)) #  T2D genes in my pool
  #   myN <- length(genome) - length(which(myGeneList %in% genome))  # genes in genome not T2D 
  #   myK <- length(myGeneSet) # size of set tested
  #   return(dhyper(x=myX, m=myM, n=myN, k=myK))
  # }
  
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
  
  for(z in stage){
    
    # p-value for enrichment in T2D genes using my differentially expressed genes for each stage
    prob_results[[z]]=1-phyper(q=length(which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID))-1,
                               m=length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id)),
                               n=length(dge_cc$genes$ensembl_gene_id) - length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id)),
                               k=length(sig_stages[[z]]$GeneID))
    
    
    # same, for random set of genes of same size. random gene set distribution of p-values
    probab_random=numeric()
    
    
    for(i in 1:10000){
      probab_random[i] <- phyperRandom( test[[d]]$GeneID, sig_stages[[z]]$GeneID, dge_cc$genes$ensembl_gene_id )  # get probs from random gene sets, for each diff expr set tested
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
  
  df_list[[d]] <- do.call("rbind", sums_probs)
  write.table(df_list[[d]],quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                                                        trait,
                                                                        "/",
                                                                        currentDate,
                                                                        "_all_DEA_timecourse_pvals_hyperg_",
                                                                        trait,
                                                                        "_genes_enrichment_10k_permutations_",
                                                                        d,
                                                                        "_kb_from_credible_regions.txt",sep=""),sep="\t")

}

distance = c("0","50","100","200","500") # list of distances in kb

# read in to save time if I just want to plot 
df_list <- list()
for (d in distance){
  df_list[[d]]=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                     trait,
                                     "/2017-05-03_all_DEA_timecourse_pvals_hyperg_",
                                     trait,
                                     "_genes_enrichment_10k_permutations_",
                                     d,
                                     "_kb_from_credible_regions.txt",sep=""),header=T)
  
  df_list[[d]]$stage=rownames(df_list[[d]])
  df_list[[d]]$distance=rep(d,8)

}

df= do.call("rbind", df_list)
colnames(df)=c("pval","lowCI","highCI","stage","distance")

df$stage=ordered(df$stage,levels = stage)
df$distance=ordered(df$distance,levels = distance)

# plot p-values and confidence intervals
# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right
diaPalette <- c("#CADAE8", "#7883BA", "#755A91", "#CC85B1", "#C15858", "#F4B8B0", "#96665A", "#6DA567")  # Diabetologia palette

ggplot(df, aes(x=stage, y=pval, colour=distance, group=distance)) + 
  geom_errorbar(aes(ymin=lowCI, ymax=highCI), colour="black", width=.1, position=pd) +
  geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
  xlab("Stages") +
  ylab("p-value after 10k permutations") +
  scale_colour_hue(name="Distance around credible region",    # Legend label, use darker colors
                   breaks=distance,
                   labels= paste(distance,"kb",sep=" "),
                   l=40) +                    # Use darker colors, lightness=40
  ggtitle(paste("Enrichment of differentially expressed genes in ",trait," genes",sep="")) +
  expand_limits(y=0) +                        # Expand y range
  # scale_y_continuous(breaks=0:1*4) +         # Set tick every x
  theme_bw() +
  # theme(legend.justification=c(1,0),
  #       legend.position=c(1,0)) +              # Position legend in bottom right
  geom_hline(yintercept=0.05,linetype="dashed",size=0.9,col="red")

ggsave(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                  trait,
                  "/2017-05-03_all_DEA_timecourse_pvals_hyperg_",
                  trait,
                  "_genes_enrichment_10k_permutations_x_axis_stages.jpg",sep=""),
       width=10,height=8,units="in",dpi=300)


# x axis distance
# -log10

df_log=df
df_log$pval=-log10(df$pval)
df_log$lowCI=-log10(df$lowCI)
df_log$highCI=-log10(df$highCI)

colnames(df_log)[4]="Stages"
levels(df_log$Stages) <- c(levels(df_log$Stages), c("EN","BLC","GT","PF"))  # change name of stages, first increase levels of factor
df_log[which(df_log$Stages=="EN6"),"Stages"] <- "EN" 
df_log[which(df_log$Stages=="EN7"),"Stages"] <- "BLC"
df_log[which(df_log$Stages=="PGT"),"Stages"] <- "GT"
df_log[which(df_log$Stages=="PFG"),"Stages"] <- "PF"

df_log <- droplevels(df_log)  # drop unused levels (old names)

df_log$Stages <- ordered(df_log$Stages, levels = c("iPSC", "DE", "GT", "PF", "PE", "EP","EN", "BLC") ) # order levels, for colours


ggplot(df_log, aes(x=distance, y=pval, colour=Stages, group=Stages)) + 
  geom_errorbar(aes(ymin=lowCI, ymax=highCI), colour="black", width=.5, size=.9, position=pd) +
  geom_line(position=pd,size=1.5) +
  scale_colour_manual(values=diaPalette) +
  geom_point(position=pd, size=4, shape=21, fill="white") + 
  xlab("Distance (kb)") +
  ylab("Permuted p-value (-log10)") +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="red") +
  theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
           panel.background = element_blank(),
           panel.border = element_rect(fill = NA, colour = "black"), 
           legend.key = element_blank(),# legend.position = c(0.5,0.5),
           axis.title.y = element_text(face="bold", angle=90, size=16, vjust=0.2),
           axis.title.x = element_text(face="bold", size=16, vjust=0),
           axis.text.x = element_text(face="bold", colour = "black", angle=90, size=16, vjust=0.2, hjust =1 ),
           axis.text.y = element_text(face="bold", colour = "black",size=16),
           axis.ticks = element_line(colour = "black"),
           axis.line = element_line(colour = "black"),
           legend.text = element_text(face="bold", colour = "black",size=14),
           legend.title = element_text(face="bold", colour = "black",size=16))

ggsave(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                  trait,
                  "/2017-05-03_all_DEA_timecourse_pvals_hyperg_",
                  trait,
                  "_genes_enrichment_10k_permutations_logs.jpg",sep="")
       width=10,height=8,units="in",dpi=300)

# normal:

ggplot(df, aes(x=distance, y=pval, colour=stage, group=stage)) + 
  geom_errorbar(aes(ymin=lowCI, ymax=highCI), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3, shape=21, fill="white") + 
  xlab("Distance (kb) around credible region") +
  ylab("p-value after 10k permutations") +
  theme_bw() +
  geom_hline(yintercept=0.05,linetype="dashed",size=0.9,col="red")
ggsave(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                  trait,
                  "/2017-05-03_all_DEA_timecourse_pvals_hyperg_",
                  trait,
                  "_genes_enrichment_10k_permutations.jpg",sep=""),
       width=10,height=8,units="in",dpi=300)

# zoom in


ggplot(df, aes(x=distance, y=pval, colour=stage, group=stage)) + 
  geom_errorbar(aes(ymin=lowCI, ymax=highCI), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=3, shape=21, fill="white") + 
  xlab("Distance (kb) around credible region") +
  ylab("p-value after 10k permutations") +
  theme_bw() +
  geom_hline(yintercept=0.05,linetype="dashed",size=0.9,col="red") +
  ylim(0,0.08)  

ggsave(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                  trait,
                  "/2017-03-06_top500_DEA_timecourse_pvals_hyperg_",
                  trait,
                  "_genes_enrichment_10k_permutations_zoom.jpg",sep=""),
       width=10,height=8,units="in",dpi=300)
  




