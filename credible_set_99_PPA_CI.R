# making credible sets for fasting glucose

currentDate <- Sys.Date() # to save date in name of output files


cumm_99 <- function(vec){
  
  # vec is a arranged (descending) vector of PPAs
  
  # returns the indeces needed to get cummulative PPA of 0.99 (including)
  
  sum = 0
  
  count = 0
  
  out.vec <- c()
  
  for (i in 1:length(vec)){
    
    val <- vec[i]
    
    sum <- sum+val
    
    if (sum<=0.99){
      
      out.vec<-append(out.vec,i)
      
      count <- count + 1
      
    }
    
  }
  
  out.vec <- append(out.vec,count+1)
  
  return(out.vec)
  
}
# Note the function returns the indices of the vec that make up the 99% credible set, 
# the length of the returned vector is the number of SNPs in the credible set

list <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credible_set_FG_march_17/FG_sets.txt",header=T)

for(l in list$list){
credible=read.table(paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credible_set_FG_march_17/",l,sep=""),header = F)

credible_99=credible[cumm_99(credible$V8),] # return all info of credible sets
# save with name that has credible region coordinates
coordinates=credible_99$V1
chr=unique(gsub("\\:[^:]*","",coordinates))
coordinates=sort(sub(".*:","",coordinates))[c(1,length(coordinates))]
write.table(credible_99,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/credible_set_FG_march_17/",currentDate,"_credible_region_",
                                   chr,"_",paste(coordinates,collapse = "_"),".txt",sep=""),sep="\t",row.names = F,col.names=T,quote = F)


}
