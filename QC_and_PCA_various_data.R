# compare different kinds of RNA-seq data
# make PCA plot
# PCA with Xie part of paper pipeline

currentDate <- Sys.Date() # to save date in name of output files

library(limma)
library(ggbiplot)
library(ggplot2)
library(reshape2)  # to modify dataframes for ggplot2
library(readr) # fast reading of large files
library(org.Hs.eg.db)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(edgeR)

# functions


# order by stages
order_by_stages=function(counts,stage.=stage){
  nc=counts[ , grepl( "iPSC" , names( counts) ) ]  #takes columns whose name matches x
  
  #Do the same for the other stages
  
  for (s in stage.[-1])  {         
    
    i= counts[ , grepl( s , names( counts ) ) ]
    nc=cbind(nc,i)
  }
  
  return(nc)
}

# filter counts

filter_counts=function(counts){
  
  #get gene information for the Ensembl IDs
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
  all_ensembl_info <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name'),
                            filters = 'ensembl_gene_id', values = rownames(counts), mart = ensembl)
  rownames(all_ensembl_info) <- all_ensembl_info$ensembl_gene_id
  
  #check the geneID exists in Ensembl
  counts <- counts[all_ensembl_info$ensembl_gene_id,]
  #load the data into edgeR
  dge <- DGEList(counts=counts,genes=all_ensembl_info)
  
  #remove genes not expressed in at least one stage with CPM > 1
  isexpr <- rowSums(cpm(dge) > 1) >= 3   
  dge <- dge[isexpr, , keep.lib.sizes=TRUE]  # I had it previously (and always) on FALSE
  
  # keep in mind that the previous function does not check if the counts are present within the same stage, so some genes that shouldn't be here remain
  # from what I've seen, they should be just a handful
  
  
  #remove genes that are not annotated as lincRNA or protein_coding
  dge <- dge[which(dge$genes$gene_biotype == "protein_coding" | dge$genes$gene_biotype == "lincRNA"),]
  #use only autosomal genes
  chr <- c(1:22)
  dge <- dge[which(dge$genes$chromosome_name %in% chr),]
  return(dge)
}

# plot SDC

plot_sdc=function(x){
  sampleDists <- dist(t(x))   
  #This function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix.
  
  
  # by default, "euclidean"
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(stages, samples,c(rep("new",24),rep("old",12)), sep=" - ")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  png("/Users/Marta/Documents/WTCHG/DPhil/Plots/Distance_clustering_new_vs_old_diff_commongenes_75bp.png",
      width=10,height=8,units="in",res=300,pointsize = 13)
  
  sdc= pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
  
  dev.off()
  
  return(sdc)
}

################## elements to be plotted

donors= c("Ad2.1","Ad3.1","Neo1.1")  # data comes from three donors
stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","ENstage6", "ENstage7") # 8 stages
old_stages=c("iPSC","DE","PGT","PFG","PE","ENstage6") # 6 stages
old_donors= c("Ad2.1","Ad3.4")
Xie_stages=c("ES","DE","PGT","PFG","PE","late_PE","polyhormonal","matured_in_vivo")
Xie_donors=c(rep("Xie1",10),rep("Xie2",5),rep("Xie3",6),rep("Xie4",4))
Nica_stages=c("beta","nonbeta","islet")
# Nica_donors=c(paste("ID",1:11,sep=""),
#               paste("ID",c(2,4,6,9,11),sep=""),
#               paste("ID",c(2,4,5,6,8,9,11),sep=""))
Nica_donors=c(rep("adult",23))
##########
# load new differentiated cells' data

cc=read_tsv("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/31.01.2017.Differentiation_v2.gene.counts.tsv") # load unfiltered conservative counts
cc=as.data.frame(cc)
rownames(cc)=cc[,1] # make gene id rownames
cc=cc[,c(-1:-2)] # take out first two columns
rownames(cc) <- gsub("(ENSG[0-9]+).*", "\\1", rownames(cc)) # cut everything after the dot
cc=order_by_stages(cc) # order 


colnames(cc)=paste( rep(stage,each=3),rep(donors,8),sep="_" )  # rename columns

# old iPSC data Islets 2016
#75 bp read
old_diff=read.table("/Users/marta/Documents/WTCHG/DPhil/Data/Diff_v2/Martijn_paper_newcounts/trim_75bp/09.03.2017.trim_75bp.gene.counts.tsv",header=T,check.names=F,row.names=1)

old_diff <- old_diff[,13:ncol(old_diff)]
#cut the version number of the Ensembl IDs
rownames(old_diff) <- gsub("(ENSG[0-9]+).*", "\\1", rownames(old_diff))

colnames(old_diff)=paste(rep(old_stages,2),rep(old_donors,each=6),sep="_")  # change column names

old_diff=order_by_stages(old_diff,stage=old_stages)


# Xie 2013

Xiecounts <- read.table("/Users/marta/Documents/WTCHG/DPhil/Data/Reference/Xie_2013_and_Nica/201215_stembancc.xie.nica.gene.counts",
                        header=T, check.names=F,row.names=1)

#cut the version number of the Ensembl IDs
rownames(Xiecounts) <- gsub("(ENSG[0-9]+).*", "\\1", rownames(Xiecounts))
#get gene information for the Ensembl IDs
ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
all_ensembl_info <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name'),
                          filters = 'ensembl_gene_id', values = rownames(Xiecounts), mart = ensembl)
rownames(all_ensembl_info) <- all_ensembl_info$ensembl_gene_id

#check the geneID still exists in Ensembl
Xiecounts <- Xiecounts[all_ensembl_info$ensembl_gene_id,]

# remove unwanted samples:

Xiecounts=Xiecounts[,c(1,2,26:73)]


Xie <- c("DP1_ES_1","DP1_DE_1","DP1_GT_1","DP1_PF_1","DP1_PE_1",
         "DP1_ES_2","DP1_DE_2","DP1_GT_2","DP1_PF_2","DP1_PE_2",
         "DP3_ES","DP3_DE","DP3_GT","DP3_PF","DP3_PE",
         "DP3-4-5_late_PE_1","DP3-4-5_late_PE_2","DP3-4-5_late_PE_3",
         "DP3-4-5_poly_1","DP3-4-5_poly_2","DP3-4-5_poly_3", # PH
         "E2147_in-vivo_matured_1","E2147_in-vivo_matured_2","E2182_in-vivo_matured_1","E2182_in-vivo_matured_2") # FE
Xiecounts=Xiecounts[,Xie]  # sort

samples <- c(rep(donors,8),rep(old_donors,6),rep(Xie_donors))   # create the design matrix
samples=as.factor(samples)
stages <-c(rep(stage,each=3),rep(old_stages,each=2),
           rep(Xie_stages[1:5],3),rep(Xie_stages[6],3),rep(Xie_stages[7],3),rep(Xie_stages[8],4))
stages <- as.factor(stages)
study <- c(rep("new",24),rep("old",12),rep("Xie",25))
design <- model.matrix(~stages + samples + study)   # I'm not sure if I should group later stages than DE together.

combined_commongenes <- merge(cc, old_diff, by=0, all=TRUE)  # merge by row names (by=0 or by="row.names")
rownames(combined_commongenes) = combined_commongenes$Row.names  # give row names again
combined_commongenes <- combined_commongenes[,c(2:ncol(combined_commongenes))] # remove first column
combined_commongenes <- merge(combined_commongenes, Xiecounts, by=0, all=TRUE)  # can only merge 2 dataframes at any time

combined_commongenes=na.omit(combined_commongenes)   # remove rows that contain NA values (reduces table to size of smallest table)
rownames(combined_commongenes)=combined_commongenes[,1] # make gene id rownames
combined_commongenes=combined_commongenes[,-1] # take out first column

filtered_combined_commongenes=filter_counts(combined_commongenes)   # 16013 genes remaining 

filtered_combined_commongenes <- calcNormFactors(filtered_combined_commongenes)    # Calculate normalization factors. TMM by default
# save as an object for later analyses
#save(filtered_combined_commongenes, file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_old_and_new_filtered_75bp.xz" , compress="xz")
# save as table for paper
x=cbind(filtered_combined_commongenes$genes,filtered_combined_commongenes$counts)
write.table(x, file=paste("/Users/marta/Documents/WTCHG/DPhil/Data/Diff_v2/",
                          Sys.Date(),"_filtered_counts_old_new_Xie.txt",sep=""),
            quote = F,row.names = F,col.names = T,sep="\t")
rm(x)
v <- voom(filtered_combined_commongenes,design,plot=F) # voom normalize the read counts


# plot PCA
plot_pca_1=function(x,s=samples,st=stages){ # color: stage # shape: samples
  pca1<-prcomp(t(x), retx=TRUE)
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,Samples=samples,Stages=stages,Experiment=c(rep("new",24),rep("old",12),rep("Xie",25)))
  levels(pcs$Stages) <- c(levels(pcs$Stages), c("EN","BLC","GT","PF","late PE","matured in vivo"))  # change name of stages, first increase levels of factor
  pcs[which(pcs$Stages=="ENstage6"),"Stages"] <- "EN" 
  pcs[which(pcs$Stages=="ENstage7"),"Stages"] <- "BLC"
  pcs[which(pcs$Stages=="PGT"),"Stages"] <- "GT"
  pcs[which(pcs$Stages=="PFG"),"Stages"] <- "PF"
  pcs[which(pcs$Stages=="late_PE"),"Stages"] <- "late PE"
  pcs[which(pcs$Stages=="matured_in_vivo"),"Stages"] <- "matured in vivo"
  
  pcs <- droplevels(pcs)  # drop unused levels (old names)
  
  pcs$Stages <- ordered(pcs$Stages, levels = c("iPSC", "ES","DE", "GT", "PF", "PE", "EP",
                                               "late PE","EN","polyhormonal", "BLC","matured in vivo") ) # order levels, for colours
  
  
  diaPalette <- c("#CADAE8","#CADAE1", "#7883BA", "#755A91", "#CC85B1", "#C15858", "#F4B8B0",
                  "#C15818","#96665A","#96165A", "#6DA567","#6DE567")  # Diabetologia palette
  
  p <- ggplot(pcs,aes(x=PC1,y=PC2, color=Stages, shape=Samples)) +
    geom_point(size=4, aes(fill=Stages), 
                           #alpha=as.character(Experiment)),
               stroke=1) +
    geom_point(size=2.5,stroke=1.5) +      
    scale_shape_manual(values=c(22,24,25,21,1,2,3,4)) +
    # scale_alpha_manual(values=c("Old"=0, "New"=diaPalette),name="Experiment",guide="none") + # guide="none" takes out legend for this alpha
    scale_colour_manual(values=diaPalette) +
    scale_fill_manual(values=diaPalette) +
    xlab (paste0( "PC1: " ,percentVar[ 1 ],"% variance")) + 
    ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  
  
  p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
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
  
  plot(p)
  # ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_PCA.tiff",p,compression="lzw")
  return(p)
}

plot_pca_2=function(x,s=samples,st=stages){ # color: experiment # shape: stage -- drop samples
  pca1<-prcomp(t(x), retx=TRUE)
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,Samples=samples,Stages=stages,Experiment=c(rep("Current",24),rep("Previous",12),rep("Xie",25)))
  levels(pcs$Stages) <- c(levels(pcs$Stages), c("EN","BLC","GT","PF","late PE","matured in vivo"))  # change name of stages, first increase levels of factor
  pcs[which(pcs$Stages=="ENstage6"),"Stages"] <- "EN" 
  pcs[which(pcs$Stages=="ENstage7"),"Stages"] <- "BLC"
  pcs[which(pcs$Stages=="PGT"),"Stages"] <- "GT"
  pcs[which(pcs$Stages=="PFG"),"Stages"] <- "PF"
  pcs[which(pcs$Stages=="late_PE"),"Stages"] <- "late PE"
  pcs[which(pcs$Stages=="matured_in_vivo"),"Stages"] <- "matured in vivo"
  
  pcs <- droplevels(pcs)  # drop unused levels (old names)
  
  pcs$Stages <- ordered(pcs$Stages, levels = c("iPSC", "ES","DE", "GT", "PF", "PE","late PE", "EP",
                                               "EN","polyhormonal", "BLC","matured in vivo") ) # order levels, for colours
  
  diaPalette <- c("#96165A","#755A91","#CADAE1")  # Diabetologia palette
  
  p <- ggplot(pcs,aes(x=PC1,y=PC2, color=Experiment, shape=Stages)) +
    geom_point(size=4, 
               stroke=1.5) +
    scale_shape_manual(values=c(16, # filled circle
                                1, # empty circle
                                8, # *
                                4, # x
                                10, # circle + 
                                9,# diamond +
                                2, # empty triangle
                                3, # +
                                15, # filled square
                                0, # empty square
                                14, # square triangle
                                12)) + # square + 
    scale_alpha_manual(values=c("Current"=diaPalette, "Previous"=diaPalette,"Xie"=0),name="Experiment",guide="none") + # guide="none" takes out legend for this alpha
    #scale_shape_manual(values=1:nlevels(pcs$Stages))  +
    scale_colour_manual(values=diaPalette) +
    xlab (paste0( "PC1: " ,percentVar[ 1 ],"% variance")) + 
    ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  
  
  p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    panel.border = element_rect(fill = NA, colour = "black"), 
                    #legend.key = element_blank(),# legend.position = c(0.5,0.5),
                    axis.title.y = element_text(face="bold", angle=90, size=16, vjust=0.2),
                    axis.title.x = element_text(face="bold", size=16, vjust=0),
                    axis.text.x = element_text(face="bold", colour = "black", angle=90, size=16, vjust=0.2, hjust =1 ),
                    axis.text.y = element_text(face="bold", colour = "black",size=16),
                    axis.ticks = element_line(colour = "black"),
                    axis.line = element_line(colour = "black"),
                    legend.text = element_text(face="bold", colour = "black",size=14),
                    legend.title = element_text(face="bold", colour = "black",size=16,vjust = -1),
                    legend.key = element_rect(size = 0.7,color = "white",fill="white"), # give more space to elements within legend
                    legend.key.size = unit(2.1, 'lines'))
  
  plot(p)
  # ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_PCA.tiff",p,compression="lzw")
  return(p)
}
plot_pca_3=function(x,s=samples,st=stages){ # color: experiment # shape: stage -- drop samples
  pca1<-prcomp(t(x), retx=TRUE)
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,Samples=samples,Stages=stages,Experiment=c(rep("Current",24),rep("Previous",12),rep("Xie",25)))
  levels(pcs$Stages) <- c(levels(pcs$Stages), c("EN","BLC","GT","PF","late PE","matured in vivo"))  # change name of stages, first increase levels of factor
  pcs[which(pcs$Stages=="ENstage6"),"Stages"] <- "EN" 
  pcs[which(pcs$Stages=="ENstage7"),"Stages"] <- "BLC"
  pcs[which(pcs$Stages=="PGT"),"Stages"] <- "GT"
  pcs[which(pcs$Stages=="PFG"),"Stages"] <- "PF"
  pcs[which(pcs$Stages=="late_PE"),"Stages"] <- "late PE"
  pcs[which(pcs$Stages=="matured_in_vivo"),"Stages"] <- "matured in vivo"
  
  pcs <- droplevels(pcs)  # drop unused levels (old names)
  
  pcs$Stages <- ordered(pcs$Stages, levels = c("iPSC", "ES","DE", "GT", "PF", "PE","late PE", "EP",
                                               "polyhormonal","EN", "BLC","matured in vivo") ) # order levels, for colours
  
  diaPalette <- c("#96165A","#755A91","#CADAE1")  # Diabetologia palette
  
  p <- ggplot(pcs,aes(x=PC1,y=PC2, color=Experiment, shape=Stages,fill=Experiment)) +
    geom_point(size=4, 
               aes(fill=Experiment),
               stroke=1.5) +
    scale_shape_manual(values=c(21, 
                                21, 
                                10,
                                22, 
                                23, 
                                24, 
                                25,
                                3, 
                                4, 
                                8, 
                                9, 
                                14)) +  
    scale_fill_manual(values=c("Current"="#96165A", "Previous"="#755A91","Xie"="#FFFFFF"),name="Experiment",guide="none") + # guide="none" takes out legend for this alpha
    scale_colour_manual(values=diaPalette) +
    xlab (paste0( "PC1: " ,percentVar[ 1 ],"% variance")) + 
    ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  
  
  p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    panel.border = element_rect(fill = NA, colour = "black"), 
                    #legend.key = element_blank(),# legend.position = c(0.5,0.5),
                    axis.title.y = element_text(face="bold", angle=90, size=16, vjust=0.2),
                    axis.title.x = element_text(face="bold", size=16, vjust=0),
                    axis.text.x = element_text(face="bold", colour = "black", angle=90, size=16, vjust=0.2, hjust =1 ),
                    axis.text.y = element_text(face="bold", colour = "black",size=16),
                    axis.ticks = element_line(colour = "black"),
                    axis.line = element_line(colour = "black"),
                    legend.text = element_text(face="bold", colour = "black",size=14),
                    legend.title = element_text(face="bold", colour = "black",size=16,vjust = -1),
                    legend.key = element_rect(size = 0.7,color = "white",fill="white"), # give more space to elements within legend
                    legend.key.size = unit(2.1, 'lines'))
  
  plot(p)
  # ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_PCA.tiff",p,compression="lzw")
  return(p)
}

plot_pca_4=function(x,s=samples,st=stages){ # color: experiment # shape: stage -- drop samples
  pca1<-prcomp(t(x), retx=TRUE)
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,Samples=samples,Stages=stages,Experiment=c(rep("Current",24),rep("Previous",12),rep("Xie",25)))
  levels(pcs$Stages) <- c(levels(pcs$Stages), c("EN","BLC","GT","PF","late PE","matured in vivo"))  # change name of stages, first increase levels of factor
  pcs[which(pcs$Stages=="ENstage6"),"Stages"] <- "EN" 
  pcs[which(pcs$Stages=="ENstage7"),"Stages"] <- "BLC"
  pcs[which(pcs$Stages=="PGT"),"Stages"] <- "GT"
  pcs[which(pcs$Stages=="PFG"),"Stages"] <- "PF"
  pcs[which(pcs$Stages=="late_PE"),"Stages"] <- "late PE"
  pcs[which(pcs$Stages=="matured_in_vivo"),"Stages"] <- "matured in vivo"
  
  pcs <- droplevels(pcs)  # drop unused levels (old names)
  
  pcs$Stages <- ordered(pcs$Stages, levels = c("iPSC", "ES","DE", "GT", "PF", "PE","late PE", "EP",
                                               "polyhormonal","EN", "BLC","matured in vivo") ) # order levels, for colours
  
  diaPalette <- c("#96165A","#755A91","#CADAE1")  # Diabetologia palette
  
  p <- ggplot(pcs,aes(x=PC1,y=PC2, color=Experiment, shape=Stages,fill=Experiment)) +
    geom_point(size=4, 
               aes(fill=Experiment),
               stroke=1.5) +
    scale_shape_manual(values=c(3, 
                                3, 
                                24,
                                22, 
                                14, 
                                25, 
                                10,
                                21, 
                                4, 
                                8, 
                                9, 
                                23)) +  
    scale_fill_manual(values=c("Current"="#96165A", "Previous"="#755A91","Xie"="#FFFFFF"),name="Experiment",guide="none") + # guide="none" takes out legend for this alpha
    scale_colour_manual(values=diaPalette) +
    xlab (paste0( "PC1: " ,percentVar[ 1 ],"% variance")) + 
    ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  
  
  p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    panel.border = element_rect(fill = NA, colour = "black"), 
                    #legend.key = element_blank(),# legend.position = c(0.5,0.5),
                    axis.title.y = element_text(face="bold", angle=90, size=16, vjust=0.2),
                    axis.title.x = element_text(face="bold", size=16, vjust=0),
                    axis.text.x = element_text(face="bold", colour = "black", angle=90, size=16, vjust=0.2, hjust =1 ),
                    axis.text.y = element_text(face="bold", colour = "black",size=16),
                    axis.ticks = element_line(colour = "black"),
                    axis.line = element_line(colour = "black"),
                    legend.text = element_text(face="bold", colour = "black",size=14),
                    legend.title = element_text(face="bold", colour = "black",size=16,vjust = -1),
                    legend.key = element_rect(size = 0.7,color = "white",fill="white"), # give more space to elements within legend
                    legend.key.size = unit(2.1, 'lines'))
  
  plot(p)
  # ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_PCA.tiff",p,compression="lzw")
  return(p)
}


p1=plot_pca_1(v$E)

p2=plot_pca_2(v$E)

ggsave(p1,"/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_diff_vs_old_75bp_Xie_commongenes_PC1and2.jpg",width=10,height=8,units="in",dpi=300)
ggsave(p2,"/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_diff_vs_old_75bp_Xie_commongenes_PC1and2_dropsample.jpg",width=10,height=8,units="in",dpi=300)


batch_corrected= removeBatchEffect(v$E,study)


p1_corrected=plot_pca_1(batch_corrected)
p2_corrected=plot_pca_2(batch_corrected)
p3_corrected=plot_pca_3(batch_corrected)
p4_corrected=plot_pca_4(batch_corrected)

ggsave(p1_corrected,"/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_diff_vs_old_75bp_Xie_commongenes_PC1and2_batchcorrected.jpg",width=10,height=8,units="in",dpi=300)
ggsave(p2_corrected,file="/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_diff_vs_old_75bp_Xie_commongenes_PC1and2_batchcorrected_dropsample.jpg",width=11,height=8,units="in",dpi=500)
ggsave(p3_corrected,file="/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_diff_vs_old_75bp_Xie_commongenes_PC1and2_batchcorrected_dropsample_Xienofill.jpg",width=11,height=8,units="in",dpi=500)
ggsave(p4_corrected,file="/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_diff_vs_old_75bp_Xie_commongenes_PC1and2_batchcorrected_dropsample_Xienofill2.jpg",width=11,height=8,units="in",dpi=500)

# Xie and Nica


Nicacounts <- read.table("/Users/marta/Documents/WTCHG/DPhil/Data/Reference/Xie_2013_and_Nica/201215_stembancc.xie.nica.gene.counts",
                        header=T, check.names=F,row.names=1)

#cut the version number of the Ensembl IDs
rownames(Nicacounts) <- gsub("(ENSG[0-9]+).*", "\\1", rownames(Nicacounts))
#get gene information for the Ensembl IDs
ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
all_ensembl_info <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name'),
                          filters = 'ensembl_gene_id', values = rownames(Nicacounts), mart = ensembl)
rownames(all_ensembl_info) <- all_ensembl_info$ensembl_gene_id

#check the geneID still exists in Ensembl
Nicacounts <- Nicacounts[all_ensembl_info$ensembl_gene_id,]

# remove unwanted samples:

Nicacounts=Nicacounts[,c(1,2,26:73)]

Nica<- c( "ID1beta","ID2beta","ID3beta","ID4beta","ID5beta", "ID6beta",
         "ID7beta", "ID8beta", "ID9beta","ID10beta","ID11beta",
         "ID2nonbeta","ID4nonbeta", "ID6nonbeta", "ID9nonbeta","ID11nonbeta",
         "ID2islet","ID4islet","ID5islet","ID6islet","ID8islet","ID9islet","ID11islet") 
Nicacounts=Nicacounts[,Nica]  # sort


samples <- c(rep(donors,8),rep(old_donors,6),Xie_donors,Nica_donors)   # create the design matrix
samples=as.factor(samples)
stages <-c(rep(stage,each=3),rep(old_stages,each=2),
           rep(Xie_stages[1:5],3),rep(Xie_stages[6],3),rep(Xie_stages[7],3),rep(Xie_stages[8],4),
           rep(Nica_stages[1],11),rep(Nica_stages[2],5),rep(Nica_stages[3],7))
stages <- as.factor(stages)
study <- c(rep("new",24),rep("old",12),rep("Xie",25),rep("Nica",23))
design <- model.matrix(~stages + samples + study) 


combined_commongenes <- merge(cc, old_diff, by=0, all=TRUE)  # merge new and old
rownames(combined_commongenes) = combined_commongenes$Row.names  # give row names again
combined_commongenes <- combined_commongenes[,c(2:ncol(combined_commongenes))] # remove first column
combined_commongenes <- merge(combined_commongenes, Xiecounts, by=0, all=TRUE)  # merging Xie
rownames(combined_commongenes) = combined_commongenes$Row.names  # give row names again
combined_commongenes <- combined_commongenes[,c(2:ncol(combined_commongenes))] # remove first column
combined_commongenes <- merge(combined_commongenes, Nicacounts, by=0, all=TRUE)  # merging Nica
combined_commongenes=na.omit(combined_commongenes)   # remove rows that contain NA values (reduces table to size of smallest table)
rownames(combined_commongenes)=combined_commongenes[,1] # make gene id rownames
combined_commongenes=combined_commongenes[,-1] # take out first column

filtered_combined_commongenes=filter_counts(combined_commongenes)   # 16504 genes remaining 

filtered_combined_commongenes <- calcNormFactors(filtered_combined_commongenes)    # Calculate normalization factors. TMM by default
# save as an object for later analyses
#save(filtered_combined_commongenes, file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_old_and_new_filtered_75bp.xz" , compress="xz")

v <- voom(filtered_combined_commongenes,design,plot=F) # voom normalize the read counts


# plot PCA
plot_pca=function(x,s=samples,st=stages){
  pca1<-prcomp(t(x), retx=TRUE)
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,Samples=samples,Stages=stages,
               Experiment=c(rep("new",24),rep("old",12),rep("Xie",25),rep("Nica",23)))
  levels(pcs$Stages) <- c(levels(pcs$Stages), c("EN","BLC","GT","PF"))  # change name of stages, first increase levels of factor
  pcs[which(pcs$Stages=="ENstage6"),"Stages"] <- "EN" 
  pcs[which(pcs$Stages=="ENstage7"),"Stages"] <- "BLC"
  pcs[which(pcs$Stages=="PGT"),"Stages"] <- "GT"
  pcs[which(pcs$Stages=="PFG"),"Stages"] <- "PF"
  
  pcs <- droplevels(pcs)  # drop unused levels (old names)
  
  pcs$Stages <- ordered(pcs$Stages, levels = c("iPSC", "ES","DE", "GT", "PF", "PE", "EP",
                                               "late_PE","EN","polyhormonal", "BLC","matured_in_vivo",
                                               "beta","nonbeta","islet") ) # order levels, for colours
  
  diaPalette <- c("#CADAE8","#CADAE1", "#7883BA", "#755A91", "#CC85B1", "#C15858", "#F4B8B0",
                  "#C15818","#96665A","#96165A", "#6DA567","#6DE567",
                  "#188914","#3576E8","#02050A")  
  
  p <- ggplot(pcs,aes(x=PC1,y=PC2, color=Stages, shape=Experiment)) +
    geom_point(size=4, aes(fill=Stages), 
               #alpha=as.character(Experiment)),
               stroke=1) +
    geom_point(size=2.5,stroke=1.5) +      
    #scale_shape_manual(values=c(22,24,25,21,1,2,3,4,5)) +
    # scale_alpha_manual(values=c("Old"=0, "New"=diaPalette),name="Experiment",guide="none") + # guide="none" takes out legend for this alpha
    scale_colour_manual(values=diaPalette) +
    scale_fill_manual(values=diaPalette) +
    xlab (paste0( "PC1: " ,percentVar[ 1 ],"% variance")) + 
    ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  
  
  p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
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
  
  plot(p)
  # ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_PCA.tiff",p,compression="lzw")
  return(p)
}

p=plot_pca(v$E)

ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_diff_vs_old_75bp_Xie_Nica_commongenes_PC2and3.jpg",p,width=10,height=8,units="in",dpi=300)


batch_corrected=removeBatchEffect(v$E,batch=study,
                                  batch2=samples,
                                  design=design[,1:15]) # remove columns we want to correct for

p_corrected=plot_pca(batch_corrected)

ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_diff_vs_old_75bp_Xie_Nica_commongenes_PC1and2_batchcorrected_batch2samples.jpg",p_corrected,width=10,height=8,units="in",dpi=300)

