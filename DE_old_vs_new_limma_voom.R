# differential expression analysis old vs new data

currentDate <- Sys.Date() # to save date in name of output files

library(limma)

load(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Diff_v2/session_objects/dge_old_and_new_filtered.xz",verbose=TRUE)  #loading the dge object for conservative counts

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



########## continue from here


design2=design[,c(1,12,2)]

v2=voom(filtered_combined_commongenes,design=design2,save.plot=TRUE) # voom normalize before fitting the linear model

fit2 <- lmFit(v2,design2)
fit2 <- eBayes(fit2)

diff_exp=topTable(fit2,coef=ncol(design2),sort.by = "none",number=nrow(fit2$coefficients))  

diff_exp_sig=diff_exp[which(diff_exp$adj.P.Val<0.01),]  #pval<0.01
