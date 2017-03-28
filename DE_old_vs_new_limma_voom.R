# differential expression analysis old vs new data

currentDate <- Sys.Date() # to save date in name of output files

library(limma)

load(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_old_and_new_filtered_75bp.xz",verbose=TRUE)  #loading the dge object for conservative counts

donors= c("Ad2.1","Ad3.1","Neo1.1")  # data comes from three donors
stage= c("iPSC", "DE", "PGT", "PFG", "PE", "EP","ENstage6", "ENstage7") # 8 stages
old_stages=c("iPSC","DE","PGT","PFG","PE","ENstage6") # 6 stages
old_donors= c("Ad2.1","Ad3.4")


# differential expression for each stage, old vs new

samples <- c(rep(donors,8),rep(old_donors,6))  
samples=as.factor(samples)
stages <-c(rep(stage,each=3),rep(old_stages,each=2))
stages <- as.factor(stages)
design <- model.matrix(~stages + samples)   # I'm not sure if I should group later stages than DE together.

experimentnew <- c(rep(1,24),rep(0,12))  # base level: old experiment
design=cbind(design,experimentnew)  # add to design


# test: old vs new, for each stage, throughout all shared stages (iPSC to EN6)
# two ways: including in the design matrix the stage vs other stages column (stage specific) -> that would give the result of old vs new for that stage
# compared to the average of all others

# then I could just take out the stages that are not the ones I'm interested in, and leave out the stage component, including just sample and experiment

# see how results change

# reorder original design matrix, and take out latest stage to make looping through column names easier:

design=cbind(design[,c(1,12,5)],c(rep(0,3),rep(1,3),rep(0,30)),design[,c(8,7,6,4,2,9,10,11)])
colnames(design)[4] <- "stagesDE"


########## 1st method
# 
# design2=design[,c(1:2,9)]   # EN6
# 
# v2=voom(filtered_combined_commongenes,design=design2,save.plot=TRUE) # voom normalize before fitting the linear model
# 
# fit2 <- lmFit(v2,design2)
# fit2 <- eBayes(fit2)
# 
# diff_exp=topTable(fit2,coef=ncol(design2),sort.by = "none",number=nrow(fit2$coefficients))  
# 
# diff_exp_sig=diff_exp[which(diff_exp$adj.P.Val<0.01),]  #pval<0.01


###### 2nd method:


diff_exp_allstages <- list()
diff_exp_allstages_nominal <- list()
for(s in old_stages){
  
  v2=voom(filtered_combined_commongenes,design=design,save.plot=TRUE) 
  rownames(design)=colnames(v2$E)   # to select rows more easily
  v3=v2[,grepl( s , colnames( v2) )]  # take only columns from selected stage
  
  design3=design[grepl( s , rownames( design) ),c(1,2)]  # same with design matrix
  
  # keep genes that have cpm>1 from voom object in all samples in either experiment old or new in the stage we are checking
  cpm_stage=cpm(filtered_combined_commongenes)  # calculate cpm
  cpm_stage=cpm_stage[,grepl(s , colnames(cpm_stage))]  # select columns from stage

  cpm_stage[cpm_stage<1]=0  # change to binary to make the selection easier
  cpm_stage[cpm_stage>=1]=1
  keep_genes=names(which(rowSums(cpm_stage[,c(1:3)])==3 | rowSums(cpm_stage[,c(4:5)])==2))  # select genes that pass filter
  v3=v3[keep_genes,]
  
  # raise filter to 10 cpm
  # # 
  # cpm_stage=cpm(filtered_combined_commongenes)  # calculate cpm
  # cpm_stage=cpm_stage[,grepl(s , colnames(cpm_stage))]  # select columns from stage
  # 
  # cpm_stage[cpm_stage<10]=0  # change to binary to make the selection easier
  # cpm_stage[cpm_stage>=10]=1
  # keep_genes=names(which(rowSums(cpm_stage[,c(1:3)])==3 | rowSums(cpm_stage[,c(4:5)])==2))  # select genes that pass filter
  # v3=v3[keep_genes,]
  # 
  ###
  fit3 <- lmFit(v3,design3)
  fit3 <- eBayes(fit3)
  
  diff_exp3=topTable(fit3,coef=ncol(design3),sort.by = "none",number=nrow(fit3$coefficients))  
  diff_exp_allstages_nominal[[s]]=diff_exp3[which(diff_exp3$P.Value<0.01 & (diff_exp3$logFC>1 | diff_exp3$logFC<(-1))),]
  diff_exp_allstages[[s]]=diff_exp3[which(diff_exp3$adj.P.Val<0.01 & (diff_exp3$logFC>1 | diff_exp3$logFC<(-1))),]  #pval<0.01 and logFC>1
   write.csv(diff_exp_allstages[[s]],file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/Diff_v2/Voom/comparison_old_and_new_study/old_with_75bp_reads/",currentDate,"_diff_exp_old_vs_new_protocol_",s,"_stage_1percent_adjpval_1cpm_old_75bp_reads.csv",sep=""),
            row.names = F)
  
}

# increasing number of significant differences between new and old expression, except for PE
# maybe due to less expressed genes?

# counting number of genes in each stage with expression >1cpm in at least two samples
zero_counts_percent = list()
for(s in old_stages){
  stage=cpm(filtered_combined_commongenes$counts)  #takes columns whose name matches x
  stage=stage[ , grepl( s , colnames( stage) ) ]
  stage[stage>=1]=1 # convert to table where any value above 1 cpm = 1
  stage[stage<1]=0 # convert to table where any value below 1 cpm = 0
  
  sum_of_each_row=rowSums(stage)
  zero_counts_percent[[s]]=100*(sum(sum_of_each_row<3)/length(sum_of_each_row)) # % of zero counts in at least 2 stages
}

# doesn't seem like PE has more 0 counts than others, but there are less differentially expressed genes than the average of other stages, 
# specially with the stage-specific method 


# how many logFC are negative? (that means that old exp has increased expression compared to new stage)
is_negative <- list()
for(s in old_stages){
  is_negative[[s]]=table(diff_exp_allstages[[s]]$logFC<0)
}


