# GO term association using GOStats
library(biomaRt)
library(GOstats)
library(qvalue)
library(org.Hs.eg.db)

currentDate <- Sys.Date() # to save date in name of output files

stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")


# just data used for the Differential expression analysis

load("/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_cc.xz")  # 15220 genes and lincRNA

sig_stages=list()
sig_stages_reduced <- list()
for(s in stage){
  # read in DEA data for each stage
  sig_stages[[s]] <- read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages/2017-03-01_sig_",s,"_diff_expression_maxvals_across-stages_results_logFC1.csv",sep=""),header = T) 
  sig_stages[[s]]=sig_stages[[s]][ order(sig_stages[[s]]$adj.P.Val,decreasing=F),] # order by adjusted p.value of "timecourse" (across-stages) method
  # top 500 genes
  sig_stages_reduced[[s]]=sig_stages[[s]][c(1:500),]
}



#perform GO enrichment using GOstats

ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")

##get the entrez IDs for the genes of interest. Takes Ensembl gene ids from my data
all_ensembl_to_entrez <- getBM(attributes=c('ensembl_gene_id', 'entrezgene'), 
                               filters = 'ensembl_gene_id', 
                               values = rownames(dge_cc), mart = ensembl)  # all genes from dge object

sig_genes <- list()

for(i in stage){
  
  sig_genes[[i]] <- getBM(attributes=c('ensembl_gene_id', 'entrezgene'), 
                          filters = 'ensembl_gene_id', 
                          values = sig_stages[[i]]$ensembl_gene_id, mart = ensembl)
}

# results: data frame with ensembl_gene_id and entrez gene id


entrez_object <- org.Hs.egGO
mapped_genes <- mappedkeys(entrez_object) # all entrez keys that have GO id
all_with_go <- unique(intersect(all_ensembl_to_entrez$entrezgene,mapped_genes)) # all entrez ids that have GO id from dge object

# now for every sig DEA results
sig_genes_with_go <- list()

for(i in stage){
  sig_genes_with_go[[i]]=unique(intersect(sig_genes[[i]]$entrezgene,mapped_genes))
}

GO_list <- list()

for (i in stage){
  params <- new('GOHyperGParams',
                geneIds=sig_genes_with_go[[i]],  # list of genes I'm testing
                universeGeneIds=all_with_go, # list of genes I used in my DE analysis (all my genes with RNA-seq data)
                ontology='BP',  # biological process
                pvalueCutoff=1,  # p value=1 to be able to calculate FDR adjusted p values later (q-values)
                conditional=T,
                testDirection='over',
                annotation="org.Hs.eg.db"
  )
  hgOver <- hyperGTest(params)
  hgOver  
  
  result <- summary(hgOver)
  
  
  result <- result[which(result$Size > 1 & result$Count > 1),] 
  # 
  # png(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/Diff_v2/GO_histograms/",currentDate,"GO_histogram_",i,".png",sep=""), type="cairo",
  #     width=10,height=10,units="in",res=300,pointsize = 13)
  # 
  # hist(result$Pvalue, 
  #      nclass = 20,
  #      main = paste("Histogram of GO terms (BP) pvals:" , i, "stage")
  # )  #inspect histogram of p-values
  # dev.off()
  # save
  
  # if x -> q value, if not p.adjust
  # U shaped. Can't use qvalue package. Use p.adjust function using BH method.
  
    # result <- cbind(result,Qvalue=qvalue(result$Pvalue)$qvalue)
    # result <- result[which(result$Qvalue<0.05),]
    # 
  # # ## calculate adjusted pvalues (BH)
   result <- cbind(result,adjust_pval=p.adjust(result$Pvalue,"BH"))
   result <- result[which(result$adjust_pval<0.05),] # saving GO terms with significant adjusted p-values

  GO_list[[i]]<- result
  
}

for(i in stage){
  # change name according to method used to calculate adjusted p-vals and logFC>1 or 0
  write.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/across-stages/GOstats/",currentDate,"_GOstats_GO_BP_",i,"_sig_result_BH_adjustment_logFC1.csv",sep=""),GO_list[[i]],row.names=F)
}

