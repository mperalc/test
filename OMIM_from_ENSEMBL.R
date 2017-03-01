# annotate with OMIM information

library(XML)
library(romim)
library(org.Hs.eg.db)

currentDate <- Sys.Date() # to save date in name of output files


stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","EN6", "EN7")
my_key <- read.table(file="/Users/Marta/Documents/WTCHG/OMIM_key.txt")
my_key <- as.character(my_key$V1)
my_key <- set_key(my_key)    # API key for OMIM queries


sig_stages=list()

# load files with DE results

for(i in stage){
  sig_stages[[i]]= read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/2017-03-01_sig_", 
                       i,"_diff_expression_maxvals_stage-unique_and_across-stages_logFC1.csv",sep=""))
}

# get OMIM id from ensembl id
cols <- c("ENSEMBL","SYMBOL", "OMIM" )

mappedGenes=list()


for(i in stage){
 mappedGenes[[i]] <- as.data.frame(select(org.Hs.eg.db,
                                     keys=as.character(sig_stages[[i]]$ensembl_gene_id),
                                     columns=cols,
                                     keytype="ENSEMBL"))
 # takes ENSEMBL ids and retrieves OMIM ids. Some are NA (don't have associated phenotypes in OMIM)
 # also, some ENSEMBL genes are NA (doesn't recognize identifiers of some lincRNA, the protein_coding genes match)
 }
 
 #to test with just one gene:
 # mappedGenes <- as.data.frame(select(org.Hs.eg.db,
 #                                     keys="ENSG00000135446",
 #                                     columns=cols,
 #                                     keytype="ENSEMBL"))


# apply romim functions to retrieve phenotypes from my list of omim ids in every result csv table:


my_list <- lapply(mappedGenes,function(x)as.integer(x[,3]))
my_list=lapply(my_list,function(x)x[!is.na(x)])         # remove Nas

my_list_omim <-list()

for(i in stage){
print(system.time(my_list_omim[[i]] <- lapply(lapply(my_list[[i]], get_omim), get_title)))
  # get_omim gets XML result for the list of OMIM id, in order.
  # get_title gives as a list the phenotype titles from OMIM ids.
}


# save OMIM object to save time in each run 
# save(my_list_omim, file = "/Users/Marta/Documents/WTCHG/DPhil/Data/OMIM/session_objects/my_list_omim.xz" , compress="xz")
 load(file="/Users/Marta/Documents/WTCHG/DPhil/Data/OMIM/session_objects/my_list_omim.xz",verbose=TRUE)  #loading the OMIM list object

df<-list()
OMIM_unique<-list()
 
for(i in stage){
  df[[i]] = data.frame(OMIM=my_list[[i]],             # merges OMIM ids and phenotypes into dataframe
                       OMIM_phen=unlist(my_list_omim[[i]]))
  
  # merge with mappedGenes by OMIM ids
  
  OMIM_unique[[i]] <- merge(x=mappedGenes[[i]], y=df[[i]], by='OMIM')
  
  # keep NAs
  
  mappedGenes[[i]]$OMIM_term=OMIM_unique[[i]][match(mappedGenes[[i]]$OMIM,OMIM_unique[[i]]$OMIM),4]
  colnames(mappedGenes[[i]])[1]="ensembl_gene_id"
}



# merge OMIM terms result into DEA table#############################################################DO
full_tables <- list()


for(i in stage){
  
  full_tables[[i]] = merge(sig_stages[[i]],mappedGenes[[i]], by="ensembl_gene_id") 
  
 # there are repeated genes, as there might be various OMIM terms associated with a certain gene 
  write.table(full_tables[[i]],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/OMIM/",
                                                                   currentDate,"_sig_",i,"_DEA_maxvals_stage-unique_and_across-stages_logFC1_OMIM.tsv",sep=""),sep="\t")
  
}




####################################################################################################
# Creating smaller tables including just genes with OMIM terms related to pancreatic defects, development defects, 
# diabetes, obesity, endocrine problems, etc
summary_terms=c("DIABETES","GLYCEMIA","OBESITY","ENDOCRINE","PANCREAS","PANCREATIC","BODY MASS INDEX",
          "DEVELOPMENT"," LIVER ","INSULIN", "SUGAR", "GLUCAGON", "ISLET")
#liver is surrounded by spaces so as not to catch other word that include it (like "oliver")

small_table =list()
small_table_aggr= list()
for(i in stage){
  keep=list()
  
  for(s in summary_terms){
    
    
    keep[[s]]=mappedGenes[[i]][grep(s, mappedGenes[[i]]$OMIM_term), ] # look into OMIM terms and take out rows when a match for summary terms is found
    
  }
  
  keep <- do.call("rbind",keep) #combine all vectors into a matrix
  keep$summary_terms=gsub("\\..*","",rownames(keep))
  # merge DEA results into this smaller table and save
  
  small_table[[i]]=sig_stages[[i]][match(keep$ensembl_gene_id,sig_stages[[i]]$ensembl_gene_id),]
  # (check that I'm not losing OMIM terms per ENSEMBL id, as they can have more than one in the "keep" table)
  small_table[[i]]=cbind(small_table[[i]],keep[,c(3,4,5)])

  
  # fix duplicated terms here:
  small_table_subset=unique(small_table[[i]][c(1:(ncol(small_table[[i]])-1))])  # eliminated repeated rows (after taking out last column)
  

  small_table_aggr[[i]]=aggregate(cbind(OMIM,as.character(OMIM_term)) ~ ensembl_gene_id , 
                                  data = small_table_subset,paste, collapse = "/") # to aggregate OMIM-related terms by ensembl_gene_id
  
  #use merge on terms from small table I'm interested in keeping
  
  small_table[[i]] <- subset(small_table[[i]], !duplicated(ensembl_gene_id)) # taking out duplicate genes before merging
  
 
 
  small_table_aggr[[i]]=merge(small_table[[i]][c(1:9)],small_table_aggr[[i]], by="ensembl_gene_id") # merge p-vals, logFCs, etc
  colnames(small_table_aggr[[i]])[11]="OMIM_terms"
  
  
  small_table_aggr[[i]] = apply(small_table_aggr[[i]], 2, gsub, patt=",", replace=" - ") # replace commas with - for easier import in excel
  
  
  write.table(small_table_aggr[[i]],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_with_and_without_contrasts/OMIM/",
                                                                 currentDate,"_sig_",i,"_DEA_maxvals_stage-unique_and_across-stages_logFC1_OMIM_interesting_terms.tsv",sep=""),sep="\t")
}



