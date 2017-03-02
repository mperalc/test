# compare old beta-cell differentiation model (pre-2016) to new one (2016-2017).
# PCA plot

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
  
  png("/Users/Marta/Documents/WTCHG/DPhil/Plots/Distance_clustering_new_vs_old_diff_commongenes.png",
    width=10,height=8,units="in",res=300,pointsize = 13)
  
  sdc= pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
  
   dev.off()
  
  return(sdc)
}

# plot PCA
plot_pca=function(x,s=samples,st=stages){
  pca1<-prcomp(t(x), retx=TRUE)
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,sample=samples,stage=stages,experiment=c(rep("new",24),rep("old",12)))
  levels(pcs$stage) <- c(levels(pcs$stage), c("EN","BLC"))  # change name of stages, first increase levels of factor
  pcs[which(pcs$stage=="ENstage6"),"stage"] <- "EN" 
  pcs[which(pcs$stage=="ENstage7"),"stage"] <- "BLC"
  pcs <- droplevels(pcs)  # drop unused levels (old names)
  
  pcs$stage <- ordered(pcs$stage, levels = c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN", "BLC") ) # order levels, for colours
  
  
  p <- ggplot(pcs,aes(x=PC1,y=PC2, color=stage, shape=sample)) +
    geom_point(size=3, aes(fill=stage, alpha=as.character(experiment)),stroke=1) +
    geom_point(size=2.5,stroke=1.5) +      
    scale_shape_manual(values=c(22,24,25,21)) +
    scale_alpha_manual(values=c("old"=0, "new"=1),name="experiment",guide="none") + # guide="none" takes out legend for this alpha
    xlab (paste0( "PC1:" ,percentVar[ 1 ],"% variance")) + 
    ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  
  
  p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    panel.border = element_rect(fill = NA, colour = "black"), 
                    legend.key = element_blank(),# legend.position = c(0.5,0.5),
                    axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
                    axis.title.x = element_text(face="bold", size=12, vjust=0),
                    axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
                    axis.text.y = element_text(face="bold", colour = "black"),
                    axis.ticks = element_line(colour = "black"),
                    axis.line = element_line(colour = "black"))
  
  plot(p)
  # ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_PCA.tiff",p,compression="lzw")
  return(p)
}

##########
# load new differentiated cells' data

cc=read_tsv("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/31.01.2017.Differentiation_v2.gene.counts.tsv") # load unfiltered conservative counts
cc=as.data.frame(cc)
rownames(cc)=cc[,1] # make gene id rownames
cc=cc[,c(-1:-2)] # take out first two columns
rownames(cc) <- gsub("(ENSG[0-9]+).*", "\\1", rownames(cc)) # cut everything after the dot

donors= c("Ad2.1","Ad3.1","Neo1.1")  # data comes from three donors
stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","ENstage6", "ENstage7") # 8 stages
old_stages=c("iPSC","DE","PGT","PFG","PE","ENstage6") # 6 stages
old_donors= c("Ad2.1","Ad3.4")

cc=order_by_stages(cc) # order 


colnames(cc)=paste( rep(stage,each=3),rep(donors,8),sep="_" )  # rename columns


# old diff data
old_diff=read.table("/Users/marta/Documents/WTCHG/DPhil/Data/Diff_v2/Martijn_paper_newcounts/31.01.2017.Islets_paper_data_2016.gene.counts.tsv",header=T,check.names=F,row.names=1)
old_diff <- old_diff[,13:ncol(old_diff)]
#cut the version number of the Ensembl IDs
rownames(old_diff) <- gsub("(ENSG[0-9]+).*", "\\1", rownames(old_diff))

colnames(old_diff)=paste(rep(old_stages,2),rep(old_donors,each=6),sep="_")  # change column names

old_diff=order_by_stages(old_diff,stage=old_stages)

# old diff and cc have the same genes

#samples <- c(rep(donors,8),rep(old_donors,6),rep("Nica",11))   # create the design matrix
samples <- c(rep(donors,8),rep(old_donors,6))  
samples=as.factor(samples)
stages <-c(rep(stage,each=3),rep(old_stages,each=2))
stages <- as.factor(stages)
study <- c(rep("new",24),rep("old",12))
design <- model.matrix(~stages + samples + study)   # I'm not sure if I should group later stages than DE together.

combined_commongenes <- merge(cc, old_diff, by=0, all=TRUE)  # merge by row names (by=0 or by="row.names")
combined_commongenes=na.omit(combined_commongenes)   # remove rows that contain NA values (does nothing because both df have same genes)
rownames(combined_commongenes)=combined_commongenes[,1] # make gene id rownames
combined_commongenes=combined_commongenes[,-1] # take out first column

filtered_combined_commongenes=filter_counts(combined_commongenes)   # 15458 genes remaining 

filtered_combined_commongenes <- calcNormFactors(filtered_combined_commongenes)    # Calculate normalization factors. TMM by default
v <- voom(filtered_combined_commongenes,design,plot=F) # voom normalize the read counts


sdc_voom=plot_sdc(v$E) 

p=plot_pca(v$E)
p
ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_diff_vs_old_commongenes_PC1and2.jpg",p,width=10,height=8,units="in",dpi=300)

# other PCs

pca1<-prcomp(t(v$E), retx=TRUE)

percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
percentVar <- round(100 * percentVar)
pcs <- as.data.frame(pca1$x)
pcs <- cbind(pcs,sample=samples,stage=stages,experiment=c(rep("new",24),rep("old",12)))
levels(pcs$stage) <- c(levels(pcs$stage), c("EN","BLC"))  # change name of stages, first increase levels of factor
pcs[which(pcs$stage=="ENstage6"),"stage"] <- "EN" 
pcs[which(pcs$stage=="ENstage7"),"stage"] <- "BLC"
pcs <- droplevels(pcs)  # drop unused levels (old names)

pcs$stage <- ordered(pcs$stage, levels = c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN", "BLC") ) # order levels, for colours


p <- ggplot(pcs,aes(x=PC2,y=PC3, color=stage, shape=sample)) +
  geom_point(size=3, aes(fill=stage, alpha=as.character(experiment)),stroke=1) +
  geom_point(size=2.5,stroke=1.5) +      
  scale_shape_manual(values=c(22,24,25,21)) +
  scale_alpha_manual(values=c("old"=0, "new"=1),name="experiment") +
  xlab (paste0( "PC2:" ,percentVar[ 2 ],"% variance")) + 
  ylab (paste0( "PC3: ",percentVar[ 3 ],"% variance" ))


p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, colour = "black"), 
                  legend.key = element_blank(),# legend.position = c(0.5,0.5),
                  axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
                  axis.title.x = element_text(face="bold", size=12, vjust=0),
                  axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
                  axis.text.y = element_text(face="bold", colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  axis.line = element_line(colour = "black"))

plot(p)

ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_diff_vs_old_commongenes_PC2and3.jpg",p,width=10,height=8,units="in",dpi=300)


pca1<-prcomp(t(v$E), retx=TRUE)

percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
percentVar <- round(100 * percentVar)
pcs <- as.data.frame(pca1$x)
pcs <- cbind(pcs,sample=samples,stage=stages,experiment=c(rep("new",24),rep("old",12)))
levels(pcs$stage) <- c(levels(pcs$stage), c("EN","BLC"))  # change name of stages, first increase levels of factor
pcs[which(pcs$stage=="ENstage6"),"stage"] <- "EN" 
pcs[which(pcs$stage=="ENstage7"),"stage"] <- "BLC"
pcs <- droplevels(pcs)  # drop unused levels (old names)

pcs$stage <- ordered(pcs$stage, levels = c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN", "BLC") ) # order levels, for colours


p <- ggplot(pcs,aes(x=PC1,y=PC3, color=stage, shape=sample)) +
  geom_point(size=3, aes(fill=stage, alpha=as.character(experiment)),stroke=1) +
  geom_point(size=2.5,stroke=1.5) +      
  scale_shape_manual(values=c(22,24,25,21)) +
  scale_alpha_manual(values=c("old"=0, "new"=1),name="experiment") +
  xlab (paste0( "PC1:" ,percentVar[ 1 ],"% variance")) + 
  ylab (paste0( "PC3: ",percentVar[ 3 ],"% variance" ))


p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, colour = "black"), 
                  legend.key = element_blank(),# legend.position = c(0.5,0.5),
                  axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
                  axis.title.x = element_text(face="bold", size=12, vjust=0),
                  axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
                  axis.text.y = element_text(face="bold", colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  axis.line = element_line(colour = "black"))

plot(p)

ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_diff_vs_old_commongenes_PC1and3.jpg",p,width=10,height=8,units="in",dpi=300)


# 
# # top 8000 expressed genes
# 
# top_genes=filtered_combined_commongenes[order(rowSums(filtered_combined_commongenes$counts),decreasing=T),] # order by decreasing counts
# top_genes=top_genes[1:8000,] # select top 8000
# top_genes=rownames(top_genes$counts) # select gene ids
# 
# # select those genes from voom object to plot
# v2=v[which(v$genes$ensembl_gene_id %in% top_genes),]
# 
# 
# p=plot_pca(v2$E)
# p
# ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/new_diff_vs_old_top8000expr_commongenes_PC1and2.jpg",p,width=10,height=8,units="in",dpi=300)

