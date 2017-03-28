# read in eQTL info from 174 Oxford islets

# best eQTLs per gene, not all are significant after correction for multiple testing

best_eqtls <- read.table(gzfile("/Users/Marta/Documents/WTCHG/DPhil/Data/eQTLs/Oxford_174_islets_Jason/best_per_gene/permutations.all.chunks.txt.gz"), header=F)

# table info: 
  # 1. Ensembl gene id
  # 2. Number of variants tested in cis for this phenotype
  # 3. MLE of the shape1 parameter of the Beta distribution
  # 4. MLE of the shape2 parameter of the Beta distribution
  # 5. Dummy [To be described later]
  # 6. ID of the best variant found for this molecular phenotypes (i.e. with the smallest p-value)
  # 7. Distance between the molecular phenotype - variant pair
  # 8. The nominal p-value of association that quantifies how significant from 0, the regression coefficient is
  # 9. The slope associated with the nominal p-value of association [only in version > v2-184]
  # 10. A first permutation p-value directly obtained from the permutations with the direct method. 
      # This is basically a corrected version of the nominal p-value that accounts for the fact that 
      # multiple variants are tested per molecular phenotype.
  # 11. A second permutation p-value obtained via beta approximation. We advice to use this one in any downstream analysis. 

# we are only intertested in V1, V6, V8, V9, V11

best_eqtls = best_eqtls[,c(1,6,7,9,11)]

colnames(best_eqtls)=c("Ensembl_ID","SNP_ID","Distance_to_TSS","Slope","Adj.p-val")


best_eqtls$Ensembl_ID <- gsub("(ENSG[0-9]+).*", "\\1", best_eqtls$Ensembl_ID) # cut everything after the dot


# look up gene names 
ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")
all_ensembl_info <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
                          filters = 'ensembl_gene_id', values = best_eqtls$Ensembl_ID, mart = ensembl)

best_eqtls$Gene_name <- all_ensembl_info[match(best_eqtls$Ensembl_ID,all_ensembl_info$ensembl_gene_id),2]

# select SNPs with significant association with gene expression

best_eqtls = best_eqtls[which(best_eqtls$`Adj.p-val`<0.05),]  # 5% expression threshold



# all nominally significant eQTLs per gene

eqtls <- read.table(gzfile("/Users/Marta/Documents/WTCHG/DPhil/Data/eQTLs/Oxford_174_islets_Jason/all_sign_0.001/nominals_p001.all.chunks.txt.gz"), header=F)

# 1. gene id
# 2. SNP id
# 3. distance to TSS
# 4. nominal p-val
# 5. slope of linear model of association

eqtls$V1 <- gsub("(ENSG[0-9]+).*", "\\1", eqtls$V1) # cut everything after the dot


# get SNPs for gene of interest
gene_id="ENSG00000148737"  # TCF7L2

gene = eqtls[which(eqtls$V1==gene_id),]  

# select SNPs from list

SNP_list=c("rs112705080","rs547143301","rs545288070","rs560412959","rs572159979","rs144966689")
SNP = eqtls[match(SNP_list,eqtls$V2),]  

