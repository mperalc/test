
#
# INITIAL EXPLORATION OF COUNT DATA
# 
# QC assessment of count data
#
# IN: all the count data files (1 for each sample and stage)
# OUT: OUT: diff expression tables, plots of QC
#
###############################################################
#           VARIABLE 
###############################################################

# Load necessary libraries

library(readr) # fast reading of large files
library(edgeR)
library(biomaRt)  # for gene information
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(limma)
library(ggbiplot)
library(ggrepel) # provides geoms for ggplot2 to repel overlapping text labels.
library(ggplot2)
library(ggfortify)




qc_plots=TRUE     # do we want to plot the QC plots? To avoid coughing of R in the for loop for the fit
long_plots=TRUE   # same




setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/other_counts")


tpm=as.data.frame(read_tsv("27.07.2016.Differentiation_v2.gene.tpm.tsv"))   # read tpm file for plotting longitudinal tpm data


# Working directory RSEM:
setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2")

my_files <- list.files("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2",pattern=".results")   # creates list of files to import, containing gene-level count data (RSEM)
# pattern chooses only files that include ".results". That way I exclude directories I have under my folder "Diff_v2".



origin= c("Ad2.1","Ad3.1","Neo1.1")  # original samples, here called origins

stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","ENstage6", "ENstage7") # check order of stages


### 1. Load all gene count files for the different samples and stages

master=lapply(my_files,read_tsv)  # master object, is a list of matrices (24 elements, 1 per sample)
# each element of the list is a matrix with 57915 rows (gene ids) and 17 columns (gene id, length, counts, TPM, etc)
names(master)=my_files            #assigned the names of the files to each matrix


# we'll work with UNnormalized counts ("expected_count" in sample matrix)

#combine all files in one matrix


####################################################################################

# We want to use expected counts for differential expression using limma+voom. 


######### check and plot total counts in each sample, to see how comparable they are #########

try=lapply(master,"[",5)   # get ("[") 5th column (expected_counts) of each element in "counts" (list of matrices) 
z=lapply(try,sum)          # From my list of matrices (each with 1 column), I sum all the values in the column for every  element of the list
                            # This gives me the total gene count in each sample
                            #DGElist function of edgeR gives me this automatically

z=as.data.frame(unlist(z))
colnames(z)="total gene count"
# 
# par(mar = c(7, 4, 2, 2) + 0.2)
# barplot(z$`total gene count`,xlab="",ylab="total gene counts from RSEM",ylim=c(0,max(z$`total gene count`)+100000))   #plot, fix x label
# # 
# text(1:24, 
#      srt = 60, adj= 1, xpd = TRUE,
#      labels =rownames(z), cex=1)
###############################################################

### 2. Combine them in a single matrix, with genes as rows and samples-stages as columns
          ## Create sample file with info of the origin and the stage of each sample
          ## Create gene file with info of the gene_id, symbol, type, length and chromosome

# create "counts" object, with rows for genes and columns for samples. 

counts=sapply(master, function(x) x[[5]]) # get ("[") 5th column (expected_counts) of each element in "counts" (list of matrices) and combine as matrix
counts=as.data.frame(counts)
g=master$Ad2.1DE.genes.results$gene_id    # take gene ids

rownames(counts)=g

rm(g)




#cut the version number of the Ensembl IDs
rownames(counts) <- gsub("(ENSG[0-9]+).*", "\\1", rownames(counts)) # cut everything after the dot
#get gene information for the Ensembl IDs
ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
all_ensembl_info <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name'),
                          filters = 'ensembl_gene_id', values = rownames(counts), mart = ensembl)
rownames(all_ensembl_info) <- all_ensembl_info$ensembl_gene_id

counts <- counts[all_ensembl_info$ensembl_gene_id,] #keep genes in counts whose geneID still exists in Ensembl
# 52662 genes left



############## filter: remove genes that don't have counts in the 3 samples of at least 1 stage ########################


counts=counts[apply(counts[,-1], MARGIN = 1, function(x) any(x > 0)), ]  # in counts, keep rows (MARGIN=1) with at least one column larger than 0 (ignore id column)
                                                                #remaining: 30184 obs. of  24 variables
# fail= ((52662-30184)/52662)*100  #42% of the original counts are all = 0

# order columns by stages, 3 by 3. This will help the following filtering step
  
  
nc=counts[ , grepl( "iPSC" , names( counts ) ) ]  #takes columns whose name matches x
  
#Do the same for the other stages

for (s in stage[-1])  {         
    
    i= counts[ , grepl( s , names( counts ) ) ]
    nc=cbind(nc,i)
  }
  
counts=nc
rm(i,nc)

#Now to the real filtering

counts=as.matrix(counts)   # to do the following, I need a matrix


# go through the table 3 by 3
# if there are more than 2 values in three columns above 0, go to next row.
# if not, go to next 3
# if gets to end of loop and not found the values, get that row. After loop, subset matrix excluding those rows

data=counts
data[data>0]=1 # convert to table where 0=0 and 1 equals any value above 0

c=c() # initiate vector of numbers, to save rows to take out of matrix

for (r in 1:nrow(data)){
  for (i in seq(1,24,3)){   # loop through stages 3 by 3 (to jump over samples of same stage)
    
    if (sum(data[r,c(i,i+1,i+2)]) > 2) {
      
      break
    }
    if (i==22){
      
      c=c(c,r)  # if it gets to the end and there are no valid values, append to vector
    }
      
    
  }
    
}

rm(data)
# % that fail

fail=100*(length(c)/nrow(counts))   # 33%
counts=as.data.frame(counts[-c,])    # counts now has 20287 genes



# create a DGEList object using the edgeR package

all_ensembl_info=all_ensembl_info[rownames(counts),] #keep in gene info file the genes that still are in counts


group <- rep(stage,each=3)    #group by stages 
  dge <- DGEList(counts=counts,group=group,genes=all_ensembl_info) 

  
  
  #########further filtering
  
  # Usually   a gene is required to have a count of 5-10 in a library to be considered expressed in that
  # library. Users should also filter with count-per-million (CPM) rather than filtering on the
  # counts directly, as the latter does not account for differences in library sizes between samples.
  
  
  keep <- rowSums(cpm(dge)>1) >= 3          #computes cpm on dge object gene counts. 
  
  # Here, a CPM of 1 corresponds to a count of 6-7 in the smallest sample. A requirement
  # for expression in three or more libraries is used as the minimum number of samples in each
  # group is three.
  
  
  
  dge <- dge[keep, , keep.lib.sizes=FALSE]     
  
  dim(dge$counts)
  # 20273 genes left after filtering
  
  
  #remove genes that are not annotated as lincRNA or protein_coding
  dge <- dge[which(dge$genes$gene_biotype == "protein_coding" | dge$genes$gene_biotype == "lincRNA"),]
  
  #15410 genes left
  
  #use only autosomal genes (not sex genes or mitocondrial genes)
  chr <- c(1:22)
  dge <- dge[which(dge$genes$chromosome_name %in% chr),]
  
  dim(dge$counts)
  
  # 14872 genes left
  
  
  ################barplot old and new library size
  
  bars_data <-as.data.frame(dge$samples[,2])
  
  
  z1=as.data.frame(z)     # remember z contains the original library size
  
  colnames(z1)="total gene count"
  
  
  
  
  rownames(bars_data)=rownames(dge$samples)
  bars_data=merge(z1, bars_data, by=0)
  rownames(bars_data)=bars_data[,1]
  bars_data=bars_data[,-1]
  colnames(bars_data)=c("old","new")
  
  bars_data=as.data.frame(t(bars_data))
  
  nb=bars_data[ , grepl( "iPSC" , names( bars_data ) ) ]  #takes columns whose name matches x
  
  #Do the same for the other stages
  
  for (s in stage[-1])  {         
    
    i= bars_data[ , grepl( s , names( bars_data ) ) ]
    nb=cbind(nb,i)
  }
  
  bars_data=nb
  rm(i,nb)
  
  stages <- c(rep(stage,each=3))
  samples <- c(rep(c("sample1","sample2","sample3"),8))
  colnames(bars_data)=paste(stages,samples,sep="-")
  
  #plot barplot


 
  
  tiff("/Users/Marta/Documents/WTCHG/DPhil/Plots/diff_v2_lib_size_allQC_cairo_10x8in_res300_tiff.tiff", type="cairo",
       width=10,height=8,units="in",res=300,pointsize = 13,compression="lzw")   # this gives an image with size 700k, with clear text
  
  labels <- paste(stages,samples,sep="-")
  mp <- barplot(as.matrix(bars_data), main="Library size before and after trimming data", axes = FALSE, axisnames = FALSE,legend = rownames(bars_data), 
                beside=TRUE,ylim=c(0,1250000), xlab="stages",font=2)
  text(mp[2,]+0.3, par("usr")[3], labels = labels, srt = 45, adj = 1.2, cex = 0.6,   #the key to place the labels is the first argument
       xpd = TRUE)
  axis(2)
  
  dev.off()
  
  rm(bars_data)

  
  
  # 
  # bars_data=as.data.frame(t(bars_data))
  # ggplot(bars_data, aes()) + geom_bar(position="dodge")
  ######################################3
  
  
#  #normalization using TMM
#  
 dge <- calcNormFactors(dge)    # TMM by default, RLE or UQ also possible
 
 

 

#  
#  ###### to normalise, I need the design matrix of the linear model I want to fit
# 
#  
 samples <- c(rep(c("sample1","sample2","sample3"),8))   # create the design matrix
 samples=as.factor(samples)
 stages <-rep(stage,each=3) 
 stages <- as.factor(stages)
 design1 <- model.matrix(~stages + samples)                    # this matrix is X in the linear model, y=Xb+e

 
 
 # If I wanted to test 1 stage of 1 sample only vs all the others, I would need another design
d=as.data.frame(design1)
 test_iPSc=as.factor(d$stagesiPSC)    # for iPSC. Here I am testing 1 stage (with all samples) vs all the other stages
 #test_iPSC_inverse=as.factor(c(0,0,0,rep(1,21)))   #to test direction of DE
 
 
 design2=model.matrix(~samples+test_iPSc)
 
 design3= model.matrix(~test_iPSc)
 
 # voom transformation uses the experiment design matrix, and produces an EList object. Using normalised data.
 #                                   #to convert the read counts to log2-cpm, with associated weights, ready for linear modelling
 
 # #The limma-voom method assumes that rows with zero or very low counts have been removed.
 
 
 #v0=voom(dge,plot=TRUE)  # no design matrix means all stages are treated as replicates
 v1=voom(dge,design=design1,plot=TRUE)  #with all samples and stages
 
 v2=voom(dge,design=design2,plot=TRUE)     # weirdness of plot reduced after filtering using Martijn's methods. Check which one is doing it, and what does it mean
 # v3=voom(dge,design=design3,plot=TRUE)  

 
 ##################################### MDS 
 
 # Plot samples on a two-dimensional scatterplot so that distances on the plot approximate 
 #the typical log2 fold changes between the samples (MDS)
 
 mds_p=plotMDS(v2$E,gene.selection = "pairwise")    # pairwise method (default)
 mds_c=plotMDS(v2$E,gene.selection = "common")      #common method 
 
 
 
 # tiff("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2_MDS_simple_only_ours.tiff")
 # plotMDS(mds_p,labels=paste(stages, samples, sep="-"))
 # dev.off()
 
 
 # A different way of doing the same thing, to produce a nicer plot.
 
 # Rearrange data for ggplot
 
 # method: pairwise
  m_p=as.data.frame(mds_p$cmdscale.out)
  m_p <- cbind(m_p,sample=samples,stage=stages)
  colnames(m_p)=c("Dimension 1", "Dimension 2", "Samples","Stages")

  # method: common
  m_c=as.data.frame(mds_c$cmdscale.out)
  m_c <- cbind(m_c,sample=samples,stage=stages)
  colnames(m_c)=c("Dimension 1", "Dimension 2", "Samples","Stages")
 
  
  # plot pairwise
 
 tiff("/Users/Marta/Documents/WTCHG/DPhil/Plots/diff_v2_MDS_pairwise_only_ours.tiff", type="cairo",
      width=10,height=10,units="in",res=300,pointsize = 13,compression="lzw")
 ggplot(m_p) +
   geom_point(aes(`Dimension 1` ,`Dimension 2`), size = 5, color = 'grey') +
   geom_label_repel(
     aes(`Dimension 1` , `Dimension 2`, fill = factor(Samples), label = Stages),
     fontface = 'bold', color = 'white',
     box.padding = unit(0.25, "lines"),
     point.padding = unit(0.5, "lines")
   ) +
   
   theme_bw(base_size=20) + 
   theme(plot.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank() )+
   theme(panel.border= element_rect())+
   
 theme(legend.position = "bottom")  +
   ggtitle("MDS plot")+
 labs(fill = "Sample")
 dev.off()
 
 # plot common
 tiff("/Users/Marta/Documents/WTCHG/DPhil/Plots/diff_v2_MDS_common_only_ours.tiff", type="cairo",
      width=10,height=10,units="in",res=300,pointsize = 13,compression="lzw")
 ggplot(m_c) +
   geom_point(aes(`Dimension 1` ,`Dimension 2`), size = 4, color = 'grey') +
   geom_label_repel(
     aes(`Dimension 1` , `Dimension 2`, fill = factor(Samples), label = Stages),
     fontface = 'bold', color = 'white',
     box.padding = unit(0.60, "lines"),
     point.padding = unit(0.5, "lines")
   ) +
  
   theme_bw(base_size=20) + 
   theme(plot.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank() )+
   theme(panel.border= element_rect())+
   
   theme(legend.position = "bottom")  +
   ggtitle("MDS plot: common method")+
   labs(fill = "Sample")
 dev.off()
 
 
 
 
 rm(m_c,m_p,mds_p,mds_c)
 ####################################################
 
 ################## get sample distance cluster  ##################
 sampleDists <- dist(t(v2$E))   # I used the design v2, but it gives the same result with any other (depends on the expression values, v$E doesn't vary with design matrix)
 sampleDistMatrix <- as.matrix(sampleDists)
 rownames(sampleDistMatrix) <- paste(stages, samples, sep="-")
 colnames(sampleDistMatrix) <- NULL
 colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
 tiff("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2_Distance_clustering.tiff", type="cairo",
      width=10,height=8,units="in",res=300,pointsize = 13,compression="lzw") 
 pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
 dev.off()
 
 
 
 # k-means clustering (DO)
 
 # autoplot(kmeans(t(v2$E), 1), data = sampleDistMatrix, label = TRUE, label.size = 3)
 
 
 ###################get sample PCA
 
 # tried for different designs: exactly the same
 
 
 pca1<-prcomp(t(v1$E), retx=TRUE)
 
 plot(pca1, type = "l") #variance vs first 10 components
 summary(pca1)          #importance of each component (important line is "proportion of variance")
 
 percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
 percentVar <- round(100 * percentVar)
 pcs <- as.data.frame(pca1$x)
 pcs <- cbind(pcs,sample=samples,stage=stages)
 pcs$stage <- ordered(pcs$stage, levels = stage)
 
 p <- ggplot(pcs, aes(PC1, PC2, colour=stage, shape=samples)) + 
   geom_point(size = 3) + xlab (paste0( "PC1:" ,percentVar[ 1 ],"% variance")) + 
   ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
 p <- p + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.background = element_blank(),
                 panel.border = element_rect(fill = NA, colour = "black"), 
                 legend.key = element_blank(),# legend.position = c(0.5,0.5),
                 axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
                 axis.title.x = element_text(face="bold", size=12, vjust=0),
                 axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
                 axis.text.y = element_text(face="bold", colour = "black"),
                 axis.ticks = element_line(colour = "black"),
                 axis.line = element_line(colour = "black"))
 ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/diff_v2_PCA_only_ours.tiff",p,compression="lzw")
 

 
############# differential expression for all stages ################
    
   
    d=as.data.frame(design1)
    
    stagesDE=as.factor(c(0,0,0,1,1,1, rep(0,18)))  # add DE to design matrix, because the function to make it ignores this (is the baseline)
    d=cbind(d,stagesDE)
    d=d[c(1,2:8,11,9,10)] # re-order to loop through it without testing intercept or samples
    
    for (i in colnames(d[2:9])) {
      
      test=as.factor(d[,which(colnames(d)==i)]) 
      design2=model.matrix(~samples+test)
      v2=voom(dge,design=design2,save.plot=TRUE) # voom normalize before fitting the linear model
      
      
      name <- gsub("stages","",i)
      
      if(qc_plots==TRUE){
        
        tiff(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2_MVtrend",name,"vs_others.tiff"), type="cairo",
             width=10,height=10,units="in",res=300,pointsize = 13,compression="lzw") 
        plot(v2$voom.xy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", 
             pch = 16, cex = 0.25)
        title(paste("voom: Mean-variance trend:",i,"vs other stages"))
        lines(v2$voom.line, col = "red")
        dev.off()
        
      }
      
      
      #   #  After this, the usual limma pipelines for differential expression can be applied, for example:
      
      fit2 <- lmFit(v2,design2)
      fit2 <- eBayes(fit2)
      diff_exp=topTable(fit2,coef=ncol(design2),sort.by = "logFC",number=nrow(fit2$coefficients))
      
      write.csv( diff_exp, file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/voom_results",name,"vs_others_differential_expression.csv",sep="_"))
      
      
      
      plot_long=diff_exp[1:20,1:2]   # keep id and names to plot longitudinal variation of these top 20 genes 
      
      plot_long=tpm[match(plot_long$ensembl_gene_id,tpm$GeneID),]
      
      #order the columns by stage
      
      nc=plot_long[ , grepl( "iPSC" , names( plot_long ) ) ]  #takes columns whose name matches x
      
      #Do the same for the other stages
      
      for (s in stage[-1])  {         
        
        i= plot_long[ , grepl( s , names( plot_long ) ) ]
        nc=cbind(nc,i)
      }
      
      plot_long=cbind(plot_long[c(1:2)],nc)
      rm(i,nc)
      
      
      # melt data for ggplot2
      
      long = melt(plot_long, measure.vars = c(3:26))
      head(long)
      
      # rename stages and samples
      stage_2= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")  #shortening EN names
      long$variable=rep(stage_2,each=60)
      
      colnames(long)[which(names(long) == "variable")] <- "stage"
      long$sample=rep(samples,each=20)
      
      long$stage <- factor(long$stage, levels=stage_2)
      long$sample=as.factor(long$sample)
      
      
      if(long_plots==TRUE){
        tiff(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2_top20genes_long_plots",name,"vs_others.tiff",sep="_"), type="cairo",
             width=25,height=10,units="in",res=300,pointsize = 13,compression="lzw") 
        
        p <- ggplot(data = long, aes(x = stage, y = value,group=sample )) +
          ggtitle("Counts per gene, stage & sample")+
          xlab ("Differentiation stages") + 
          ylab ("Counts (TPM)") +
          geom_line(aes(col = sample), size = 1) + facet_wrap(~GeneName,scales="free")
        print(p)
        dev.off()
      }
      
      if(qc_plots==TRUE){
        
        # Volcano plot
        tiff(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2_volcano_",name,"vs_others.tiff"), type="cairo",
             width=10,height=10,units="in",res=300,pointsize = 13,compression="lzw") 
        volcanoplot(fit2,coef=ncol(design2),highlight=30,names=fit2$genes$external_gene_name,main=paste("Volcano plot(",name,"vs other stages)"))
        dev.off()
        
        # Mean-difference plot
        tiff(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2_MD_",name,"vs_others.tiff"), type="cairo",
             width=10,height=10,units="in",res=300,pointsize = 13,compression="lzw") 
        plotMD(fit2,coef=ncol(design2),status=fit2$p.value[,4]<0.001,main="Mean-difference plot",ylab = paste("log-fold change (",name,"vs other stages)"))
        dev.off()
        
        # Q-Q plot of moderated t-statistics
        tiff(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2_QQplot_mod_t_stats",name,"vs_others.tiff"), type="cairo",
             width=10,height=8,units="in",res=300,pointsize = 13,compression="lzw") 
        qqt(fit2$t[,4],df=fit2$df.residual+fit2$df.prior,main=paste("Q-Q plot of moderated t-statistics(",name,"vs other stages)"))  # careful to select the roght column in data
        abline(0,1)
        dev.off()
        
      }
      
      
      
    }   
    
    
   
   
   
   