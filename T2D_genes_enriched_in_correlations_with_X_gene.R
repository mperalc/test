# for islet diff paper 2017

# are genes/TFs within T2D loci enriched for correlation with CREBBP?

#### load in data
currentDate <- Sys.Date()

load("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_cc.xz")  # 15221 genes and lincRNA

T2D_0kb=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_T2D/2017-02-22annotated_genes_in_credible_regions_plusminus_0_kb.txt",
                   header=T)
colnames(test[[d]])[1]=c("GeneID")
test[[d]]$GeneID=as.character(test[[d]]$GeneID)


# are genes/TFs within T2D loci enriched for STRING connectivity with CREBBP?
