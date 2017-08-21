## Hypergeometric test for datasets: are they enriched in T2D/FG genes?

# in paper 2017
# for each gene Diff expr in each stage
# for each gene diff expr in the first 7 stages (below)
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

trait= "T2D"    # T2D or fasting glucose (FG)?
type="all"     #or top623
order="p-value"   # or logFC
### files

#### background "genome" data
# just data used for the Differential expression analysis

load("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_cc.xz")  # 15221 genes and lincRNA

# gene test data

distance = c("0","50","100","200","500") # list of distances in kb

test <- list()
                  

for(d in distance){
  if (trait=="FG"){
    test[[d]]=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                    trait,
                                    #"/2017-02-22annotated_genes_in_credible_regions_plusminus_",
                                    "/2017-05-03FG_annotated_genes_in_credible_regions_plusminus_",
                                    d,"_kb.txt",sep=""),header=T)
  }
  else{
    test[[d]]=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                    trait,
                                    "/2017-02-22annotated_genes_in_credible_regions_plusminus_",
                                    d,"_kb.txt",sep=""),header=T)
  
  }
  colnames(test[[d]])[1]=c("GeneID")
  test[[d]]$GeneID=as.character(test[[d]]$GeneID)
}


# files with DE results (across-stages)

# # all genes, get list of ensembl gene ids

sig_stages=list()
for(s in stage){
 if(order=="p-value"){
   sig_stages[[s]] <- read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages/2017-07-03_sig_",s,"_diff_expression_maxvals_across-stages_results_logFC1.csv",sep=""),header = T) # read in DEA data for each stage
   colnames(sig_stages[[s]])[1]="GeneID"
   sig_stages[[s]]=sig_stages[[s]][ order(sig_stages[[s]]$adj.P.Val,decreasing=F),] # order by adjusted p.value 
 }
  else{
    sig_stages[[s]] <- read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages/2017-07-03_sig_",s,"_diff_expression_maxvals_across-stages_results_logFC1.csv",sep=""),header = T) # read in DEA data for each stage
    colnames(sig_stages[[s]])[1]="GeneID"
    sig_stages[[s]]=sig_stages[[s]][ order(sig_stages[[s]]$logFC,decreasing=T),] # order by logFC 
  }
  if(type=="top623"){
    
    # top 623 genes (min significant DE genes, in PGT stage)
    sig_stages[[s]]=sig_stages[[s]][c(1:623),]
  }

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


# comparison of probabilities

phyperRandom <- function(myGeneList, myGeneSet, genome){
  myRandomGS <- sample( genome,size=length(myGeneSet) )  # take random sample of genes from whole genome, same size as set tested
  myX <- length(which(myGeneList %in% myRandomGS))  # number of T2D in random set
  myM <- length(which(myGeneList %in% genome)) #  T2D genes in my pool
  myN <- length(genome) - length(which(myGeneList %in% genome))  # genes in genome not T2D 
  myK <- length(myGeneSet) # size of set tested
  return(list(prob=1-phyper(q=myX-1,
                  m=myM, n=myN, k=myK), n.overlap=myX))
}

df_list <- list()
distance = c("0","50","100","200","500") # list of distances in kb
overlap_list <- list()  # to write down number of T2D genes in diff expr results

for(d in distance){
  
  prob_results=list()
  sums_probs=list()
  b_test <- list()
  overlap <- list()
  
  for(z in stage){
    
    # p-value for enrichment in T2D genes using my differentially expressed genes for each stage
    prob_results[[z]]=1-phyper(q=length(which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID))-1,
                               m=length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id)),
                               n=length(dge_cc$genes$ensembl_gene_id) - length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id)),
                               k=length(sig_stages[[z]]$GeneID))
    gene_list=test[[d]][which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID),]
    # # write.table(gene_list,
    #             file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
    #                   trait,
    #                   "/",
    #                   currentDate,
    #                   "_gene_list_",
    #                   trait,
    #                   "_in_top500_genes_stage_",
    #                   z,"_",d,
    #                   "_kb_around_credible_regions.txt",sep=""),
    #             sep="\t",col.names = T,row.names = F,quote = F)
    # getting amount of T2D genes in each stage and in background pool
    overlap[["total_in_background"]] = length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id))
    
    overlap[[z]] = length(which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID))
    
    # same, for random set of genes of same size. random gene set distribution of p-values
    random=list()
    random_overlaps=numeric()
    probab_random=numeric()
    
    for(i in 1:10000){
      random[[i]] <- as.data.frame(phyperRandom( test[[d]]$GeneID, sig_stages[[z]]$GeneID, dge_cc$genes$ensembl_gene_id))  # get probs from random gene sets, for each diff expr set tested
      probab_random[i] <- random[[i]]$prob
      random_overlaps[i] <- random[[i]]$n.overlap
      
    }
  
     
    # png(paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
    #                trait,
    #                "/",
    #                currentDate,
    #                "_random_overlaps_",
    #                trait,
    #                "_top500_genes_stage_",
    #                z,"_",d,
    #                "_kb_around_credible_regions.png",sep=""),
    #     type="cairo",
    #     width=8,height=5,units="in",res=300,pointsize = 12)
    hist(random_overlaps)
    # dev.off()
    #hist(probab_random)
    
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
  overlap_list[[d]] <- do.call("rbind", overlap)
  
  df_list[[d]] <- do.call("rbind", sums_probs)
  
  # ######### save tables
  # ordered by p-value
  
  # 
  # write.table(df_list[[d]],quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
  #                                                                       trait,
  #                                                                       "/",
  #                                                                       currentDate,
  #                                                                       "_top623_DEA_timecourse_pvals_hyperg_",
  #                                                                       trait,
  #                                                                       "_genes_enrichment_10k_permutations_",
  #                                                                       d,
  #                                                                       "_kb_from_credible_regions.txt",sep=""),sep="\t")
  # write.table(overlap_list[[d]],quote=F,row.names=T,col.names = F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
  #                                                                       trait,
  #                                                                       "/",
  #                                                                       currentDate,
  #                                                                       "_top623_DEA_timecourse_overlaps_hyperg_",
  #                                                                       trait,
  #                                                                       "_genes_enrichment_",
  #                                                                       d,
  #                                                                       "_kb_from_credible_regions.txt",sep=""),sep="\t")

  # ordered by logFC
  
   write.table(df_list[[d]],quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                                                        trait,
                                                                        "/",
                                                                        currentDate,
                                                                        "_top623_orderedByLogFC_DEA_timecourse_pvals_hyperg_",
                                                                        trait,
                                                                        "_genes_enrichment_10k_permutations_",
                                                                        d,
                                                                        "_kb_from_credible_regions.txt",sep=""),sep="\t")
  write.table(overlap_list[[d]],quote=F,row.names=T,col.names = F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                                                        trait,
                                                                        "/",
                                                                        currentDate,
                                                                        "_top623_orderedByLogFC_DEA_timecourse_overlaps_hyperg_",
                                                                        trait,
                                                                        "_genes_enrichment_",
                                                                        d,
                                                                        "_kb_from_credible_regions.txt",sep=""),sep="\t")
}

# 
# 
# 
# # plots ######
distance = c("0","50","100","200","500") # list of distances in kb

# read in to save time if I just want to plot
df_list <- list()
for (d in distance){
  df_list[[d]]=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                     trait,
                                     "/2017-06-01_top623_orderedByLogFC_DEA_timecourse_pvals_hyperg_",
                                   # "/2017-06-01_top623_DEA_timecourse_pvals_hyperg_",
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
# 
# ggplot(df, aes(x=stage, y=pval, colour=distance, group=distance)) +
#   geom_errorbar(aes(ymin=lowCI, ymax=highCI), colour="black", width=.1, position=pd) +
#   geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
#   xlab("Stages") +
#   ylab("p-value after 10k permutations") +
#   scale_colour_hue(name="Distance around credible region",    # Legend label, use darker colors
#                    breaks=distance,
#                    labels= paste(distance,"kb",sep=" "),
#                    l=40) +                    # Use darker colors, lightness=40
#   ggtitle(paste("Enrichment of differentially expressed genes in ",trait," genes",sep="")) +
#   expand_limits(y=0) +                        # Expand y range
#   # scale_y_continuous(breaks=0:1*4) +         # Set tick every x
#   theme_bw() +
#   # theme(legend.justification=c(1,0),
#   #       legend.position=c(1,0)) +              # Position legend in bottom right
#   geom_hline(yintercept=0.05,linetype="dashed",size=0.9,col="red")
# 
# ggsave(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
#                   trait,
#                   "/2017-05-03_top_500_DEA_timecourse_pvals_hyperg_",
#                   trait,
#                   "_genes_enrichment_10k_permutations_x_axis_stages.jpg",sep=""),
#        width=10,height=8,units="in",dpi=300)


# # x axis distance
# # -log10

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
           legend.key = element_blank(), # legend.position = c(0.5,0.5),
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
                  "/2017-06-01_top623_orderedByLogFC_DEA_timecourse_pvals_hyperg_",
                  #"/2017-06-01_top623_DEA_timecourse_pvals_hyperg_",
                  trait,
                  "_genes_enrichment_10k_permutations_logs.jpg",sep=""),
       width=10,height=8,units="in",dpi=300)

# # normal:
# 
# ggplot(df, aes(x=distance, y=pval, colour=stage, group=stage)) + 
#   geom_errorbar(aes(ymin=lowCI, ymax=highCI), colour="black", width=.1, position=pd) +
#   geom_line(position=pd) +
#   geom_point(position=pd, size=3, shape=21, fill="white") + 
#   xlab("Distance (kb) around credible region") +
#   ylab("p-value after 10k permutations") +
#   theme_bw() +
#   geom_hline(yintercept=0.05,linetype="dashed",size=0.9,col="red")
# ggsave(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
#                   trait,
#                   "/2017-05-03_top_500_DEA_timecourse_pvals_hyperg_",
#                   trait,
#                   "_genes_enrichment_10k_permutations.jpg",sep=""),
#        width=10,height=8,units="in",dpi=300)
# 
# # zoom in
# 
# 
# ggplot(df, aes(x=distance, y=pval, colour=stage, group=stage)) + 
#   geom_errorbar(aes(ymin=lowCI, ymax=highCI), colour="black", width=.1, position=pd) +
#   geom_line(position=pd) +
#   geom_point(position=pd, size=3, shape=21, fill="white") + 
#   xlab("Distance (kb) around credible region") +
#   ylab("p-value after 10k permutations") +
#   theme_bw() +
#   geom_hline(yintercept=0.05,linetype="dashed",size=0.9,col="red") +
#   ylim(0,0.08)  
# 
# ggsave(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
#                   trait,
#                   "/2017-03-06_top500_DEA_timecourse_pvals_hyperg_",
#                   trait,
#                   "_genes_enrichment_10k_permutations_zoom.jpg",sep=""),
#        width=10,height=8,units="in",dpi=300)
#   
# 
# 
# 
# 

# getting tables of genes in background


distance = c("0","50","100","200","500") # list of distances in kb

test <- list()


for(d in distance){ # FG
  test[[d]]=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_FG",
                                  "/2017-05-03FG_annotated_genes_in_credible_regions_plusminus_",
                                  d,"_kb.txt",sep=""),header=T)
  colnames(test[[d]])[1]=c("GeneID")
  test[[d]]$GeneID=as.character(test[[d]]$GeneID)
  overlap = test[[d]][which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id),]
  write.table(overlap,quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                                                        "FG",
                                                                        "/",
                                                                        currentDate,
                                                                        "_FG_genes_in background_",
                                                                        d,
                                                                        "_kb_from_credible_regions.txt",sep=""),sep="\t")
}


for(d in distance){ # T2D
  test[[d]]=read.table(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_T2D",
                                  "/2017-02-22annotated_genes_in_credible_regions_plusminus_",
                                  d,"_kb.txt",sep=""),header=T)
  colnames(test[[d]])[1]=c("GeneID")
  test[[d]]$GeneID=as.character(test[[d]]$GeneID)
  overlap = test[[d]][which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id),]
  write.table(overlap,quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                                                        "T2D",
                                                                        "/",
                                                                        currentDate,
                                                                        "_T2D_genes_in background_",
                                                                        d,
                                                                        "_kb_from_credible_regions.txt",sep=""),sep="\t")
}



# plot just significant results

# T2D

T2D=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_T2D/just_BLC_T2D.txt",
                      header=T)

T2D$pval=-log10(T2D$pval)
T2D$conf.int.low=-log10(T2D$conf.int.low)
T2D$conf.int.up=-log10(T2D$conf.int.up)
T2D$distance=as.factor(T2D$distance)
T2D$distance=ordered(T2D$distance, levels = c("0", "50", "100", "200", "500") ) # order levels
colnames(T2D)=c("distance","pval","conf.int.low","conf.int.up","Differentially expressed genes")
levels(T2D$`Differentially expressed genes`) <- c(levels(T2D$`Differentially expressed genes`), c("All","Beta cell function","Top 623"))  #  first increase levels of factor

T2D[which(T2D$`Differentially expressed genes`=="all"),"Differentially expressed genes"] <- "All"
T2D[which(T2D$`Differentially expressed genes`=="phys"),"Differentially expressed genes"] <- "Beta cell function"
T2D[which(T2D$`Differentially expressed genes`=="top_623"),"Differentially expressed genes"] <- "Top 623"
diaPalette <- c( "#7883BA",  "#C15858",  "#6DA567")  # Diabetologia palette


# ggplot(T2D, aes(x=distance, y=pval, colour=`Differentially expressed genes`)) +
#   geom_crossbar(aes(ymin=conf.int.low, ymax=conf.int.up), width=0.2,fatten = 4) +
#   #geom_line(position=pd,size=1.5) +
#   scale_colour_manual(values=diaPalette) +
#   #geom_point(position=pd, size=4, shape=21, fill="white") +
#   xlab("Distance (kB)") +
#   ylab("Permuted p-value (-log10)") +
#   theme_bw() +
#   geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="red") +
#   theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#            panel.background = element_blank(),
#            panel.border = element_rect(fill = NA, colour = "black"),
#            legend.key = element_blank(), # legend.position = c(0.5,0.5),
#            axis.title.y = element_text(face="bold", angle=90, size=16, vjust=0.2),
#            axis.title.x = element_text(face="bold", size=16, vjust=0),
#            axis.text.x = element_text(face="bold", colour = "black", angle=90, size=16, vjust=0.2, hjust =1 ),
#            axis.text.y = element_text(face="bold", colour = "black",size=16),
#            axis.ticks = element_line(colour = "black"),
#            axis.line = element_line(colour = "black"),
#            legend.text = element_text(face="bold", colour = "black",size=14),
#            legend.title = element_text(face="bold", colour = "black",size=16))
# 
# ggsave(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_T2D/just_BLC_T2D_geom_crossbar.jpg",
#        width=10,height=8,units="in",dpi=300)

ggplot(T2D, aes(x=distance, y=pval, colour=`Differentially expressed genes`,group=`Differentially expressed genes`)) +
  geom_errorbar(aes(ymin=conf.int.low, ymax=conf.int.up), width=0.1,size=1) +
  geom_line(aes(linetype = `Differentially expressed genes`),size=0.5)+
  scale_colour_manual(values=diaPalette) +
  geom_point(size = 4) +
  geom_point(size = 3, color = "white") + 
  xlab("Distance (kB)") +
  ylab("Permuted p-value (-log10)") +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="red") +
  theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           panel.border = element_rect(fill = NA, colour = "black"),
           legend.key = element_blank(), # legend.position = c(0.5,0.5),
           axis.title.y = element_text(face="bold", angle=90, size=16, vjust=0.2),
           axis.title.x = element_text(face="bold", size=16, vjust=0),
           axis.text.x = element_text(face="bold", colour = "black", angle=90, size=16, vjust=0.2, hjust =1 ),
           axis.text.y = element_text(face="bold", colour = "black",size=16),
           axis.ticks = element_line(colour = "black"),
           axis.line = element_line(colour = "black"),
           legend.text = element_text(face="bold", colour = "black",size=14),
           legend.title = element_text(face="bold", colour = "black",size=16))

ggsave(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_T2D/just_BLC_T2D_geom_point_line.jpg",
       width=12,height=8,units="in",dpi=300)

FG=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_FG/all_diff_expr_genes/FG_all_diffexpr_genes_just_BLC_EP.txt",
               header=T)

FG$pval=-log10(FG$pval)
FG$conf.int.low=-log10(FG$conf.int.low)
FG$conf.int.up=-log10(FG$conf.int.up)
FG$distance=as.factor(FG$distance)
FG$distance=ordered(FG$distance, levels = c("0", "50", "100", "200", "500") ) # order levels
colnames(FG)=c("Stage","pval","conf.int.low","conf.int.up","distance")

diaPalette <- c( "#C15858",  "#7883BA")  # Diabetologia palette

ggplot(FG, aes(x=distance, y=pval, colour=Stage,group=Stage)) +
  geom_errorbar(aes(ymin=conf.int.low, ymax=conf.int.up), width=0.1,size=1) +
  geom_line(aes(linetype = Stage),size=0.5)+
  scale_colour_manual(values=diaPalette) +
  geom_point(size = 4) +
  geom_point(size = 3, color = "white") + 
  xlab("Distance (kB)") +
  ylab("Permuted p-value (-log10)") +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05),linetype="dashed",size=0.8,col="red") +
  theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           panel.border = element_rect(fill = NA, colour = "black"),
           legend.key = element_blank(), # legend.position = c(0.5,0.5),
           axis.title.y = element_text(face="bold", angle=90, size=16, vjust=0.2),
           axis.title.x = element_text(face="bold", size=16, vjust=0),
           axis.text.x = element_text(face="bold", colour = "black", angle=90, size=16, vjust=0.2, hjust =1 ),
           axis.text.y = element_text(face="bold", colour = "black",size=16),
           axis.ticks = element_line(colour = "black"),
           axis.line = element_line(colour = "black"),
           legend.text = element_text(face="bold", colour = "black",size=14),
           legend.title = element_text(face="bold", colour = "black",size=16))

ggsave(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_FG/all_diff_expr_genes/FG_all_diffexpr_genes_just_BLC_EP.jpg",
       width=10,height=8,units="in",dpi=300)



### test for enrichment in T2D genes for each diff expr gene of the first 7 stages (all developmental)
# less strict than for each stage, as above considers mature the last stage (BLC), and the developmental
# signals are spread across 7 other stages.


df_list <- list()
distance = c("0","50","100","200","500") # list of distances in kb
overlap_list <- list()  # to write down number of T2D genes in diff expr results

# merge developmental stages in one big stage

dev=do.call("rbind", sig_stages[c(1:7)])  # dev has 7612 genes
adult=sig_stages$EN7 # adult has 1796 genes
stage=c("dev","adult")
sig_stages$dev=dev   # adding elements to list
sig_stages$adult=adult
for(d in distance){
  
  prob_results=list()
  sums_probs=list()
  b_test <- list()
  overlap <- list()
  
  for(z in stage){
    
    # p-value for enrichment in T2D genes using my differentially expressed genes for each stage
    prob_results[[z]]=1-phyper(q=length(which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID))-1,
                               m=length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id)),
                               n=length(dge_cc$genes$ensembl_gene_id) - length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id)),
                               k=length(sig_stages[[z]]$GeneID))
    gene_list=test[[d]][which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID),]
    # getting amount of T2D genes in each stage 
    write.table(gene_list,
                file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                      trait,
                      "/",
                      currentDate,
                      "_gene_list_",
                      trait,
                      "_in_",
                      type,
                      "genes_stage_",
                      z,"_",d,
                      "_kb_around_credible_regions.txt",sep=""),
                sep="\t",col.names = T,row.names = F,quote = F)
    #  and in background pool
    overlap[["total_in_background"]] = length(which(test[[d]]$GeneID %in% dge_cc$genes$ensembl_gene_id))
    
    overlap[[z]] = length(which(test[[d]]$GeneID %in% sig_stages[[z]]$GeneID))
    
    # same, for random set of genes of same size. random gene set distribution of p-values
    random=list()
    random_overlaps=numeric()
    probab_random=numeric()
    
    for(i in 1:10000){
      random[[i]] <- as.data.frame(phyperRandom( test[[d]]$GeneID, sig_stages[[z]]$GeneID, dge_cc$genes$ensembl_gene_id))  # get probs from random gene sets, for each diff expr set tested
      probab_random[i] <- random[[i]]$prob
      random_overlaps[i] <- random[[i]]$n.overlap
      
    }
    
    
    
    h=hist(random_overlaps)
    
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
  overlap_list[[d]] <- do.call("rbind", overlap)
  
  df_list[[d]] <- do.call("rbind", sums_probs)
  
  # ######### save tables
  # ordered by p-value
  

    write.table(df_list[[d]],quote=F,row.names=T,col.names = T,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                                                          trait,
                                                                          "/",
                                                                          currentDate,
                                                                          "_",
                                                                          type,
                                                                          "_ordered_by_",
                                                                          order,
                                                                          "_DEA_timecourse_pvals_hyperg_",
                                                                          trait,
                                                                          "_genes_enrichment_10k_permutations_",
                                                                          d,
                                                                          "_kb_from_credible_regions.txt",sep=""),sep="\t")
    write.table(overlap_list[[d]],quote=F,row.names=T,col.names = F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_",
                                                                          trait,
                                                                          "/",
                                                                          currentDate,
                                                                          "_",
                                                                          type,
                                                                          "_ordered_by_",
                                                                          order,
                                                                          "_DEA_timecourse_overlaps_hyperg_",
                                                                          trait,
                                                                          "_genes_enrichment_",
                                                                          d,
                                                                          "_kb_from_credible_regions.txt",sep=""),sep="\t")

}


