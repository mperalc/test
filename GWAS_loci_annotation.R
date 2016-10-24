# annotating genes around (150kb) GWAS loci
# diabetes

# as Martijn does it


############read in GWAS loci

#add and substract n kb to get chr:xxx-n:xxx+n format

## save new files

### use those files to do query


## getBM solution:

library("biomaRt")
listMarts(host="www.ensembl.org")
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")
filters = listFilters(ensembl)
filterlist<-as.list(HMapSamp$MergeCol)


results=getBM(attributes = c("hgnc_symbol","entrezgene", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),values = list(chromosomal_region=filterlist,biotype="protein_coding"), mart = ensembl)



## annotation solution:

#  use the Annotation packages rather than biomaRt, so you don't have to send a query through the wires each time. 
#First of all, install Homo sapiens from Bioconductor:
source("https://bioconductor.org/biocLite.R")
biocLite("Homo.sapiens")
library(Homo.sapiens)

# This will also install a TxDb object for homo sapiens, called TxDb.Hsapiens.UCSC.hg19.knownGene.

# To get all the gene coordinates, use the genes() function:
genes(TxDb.Hsapiens.UCSC.hg19.knownGene)


# To get the intersection with your coordinates, you must first convert them to a GenomicRanges list. 
# The easiest way is to convert the list to a dataframe, and then use the makeGRangesFromDataFrame function, 
# which should have already been loaded when with Homo.sapiens. Remember to append the prefix 'chr' to your chromosome names. Coordinates should be 1-based.
library(dplyr)
mycoords.gr = lapply(mycoords.list, function (x) {res=strsplit(x, ':')}) %>%
  unlist %>%
  as.numeric %>%
  matrix(ncol=3, byrow=T) %>%
  as.data.frame %>%
  select(chrom=V1, start=V2, end=V3) %>%
  mutate(chrom=paste0('chr', chrom)) %>%
  makeGRangesFromDataFrame


subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), mycoords.gr)



####################### Martijn's method
##look for T2D associated gene enrichment with Fisher
t2d_close <- scan("T2D_genes_closest",what="character")
t2d_50k <- scan("T2D_genes_50kb",what="character")
t2d_100k <- scan("T2D_genes_100kb",what="character")
t2d_150k <- scan("T2D_genes_150kb",what="character")

#iPSC
iPSC_t2d_close <- fisher.test(data.frame(found=c(length(which(rownames(sig_iPSC_stage) %in% t2d_close)),length(which(rownames(dge) %in% t2d_close))),
                                         notfound=c(dim(sig_iPSC_stage)[1] - length(which(rownames(sig_iPSC_stage) %in% t2d_close)),dim(dge)[1] - length(which(rownames(dge) %in% t2d_close)))))
iPSC_t2d_50 <- fisher.test(data.frame(found=c(length(which(rownames(sig_iPSC_stage) %in% t2d_50k)),length(which(rownames(dge) %in% t2d_50k))),
                                      notfound=c(dim(sig_iPSC_stage)[1] - length(which(rownames(sig_iPSC_stage) %in% t2d_50k)),dim(dge)[1] - length(which(rownames(dge) %in% t2d_50k)))))
iPSC_t2d_100 <- fisher.test(data.frame(found=c(length(which(rownames(sig_iPSC_stage) %in% t2d_100k)),length(which(rownames(dge) %in% t2d_100k))),
                                       notfound=c(dim(sig_iPSC_stage)[1] - length(which(rownames(sig_iPSC_stage) %in% t2d_100k)),dim(dge)[1] - length(which(rownames(dge) %in% t2d_100k)))))
iPSC_t2d_150 <- fisher.test(data.frame(found=c(length(which(rownames(sig_iPSC_stage) %in% t2d_150k)),length(which(rownames(dge) %in% t2d_150k))),
                                       notfound=c(dim(sig_iPSC_stage)[1] - length(which(rownames(sig_iPSC_stage) %in% t2d_150k)),dim(dge)[1] - length(which(rownames(dge) %in% t2d_150k)))))
