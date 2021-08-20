PTDO80 <- function(observed, expected, result){
  
  observed_sum <- apply(observed, 2, sum)
  max_observed <- max(observed_sum)
  observedlam <- observed_sum/max_observed
  
  ## true 0 positions figured out by MCMC
  PTDOfun <- function(thresh) {
    posi <- which(result$SZ>=thresh,TRUE)
    posirows <- NULL
    for (i in 1:(dim(posi)[1])) {
      posirows <- c(posirows,matrow(posi[i,1],posi[i,2]))
    }
    IMP <- result$Impute_All
    IMP[posirows,]<-0
    
    PTSZ=list()
    PTDO=list()
    
    for(j in 1:ncol(observed)){
      indexobserved0=(observed[,j]==0)
      predictv=IMP[,j][indexobserved0]
      Truevalue=expected[,j][indexobserved0]
      
      PTSZ[[j]]=sum(Truevalue==0 & predictv==0)/sum(Truevalue==0)
      
      PTDO[[j]]=sum(Truevalue>0 & predictv>0)/sum(Truevalue>0)
    }
    PTSZ=unlist(PTSZ)
    PTDO=unlist(PTDO)
    return(mean(PTDO))
  }
  
  threshold <- seq(0.001, max(result$SZ), length.out = 200)
  PTDOout <- unlist(lapply(threshold,PTDOfun))
  tt <- threshold[which.min(abs(PTDOout-0.8))]
  
  posi <- which(result$SZ>=tt,TRUE)
  posirows <- NULL
  for (i in 1:(dim(posi)[1])) {
    posirows <- c(posirows,matrow(posi[i,1],posi[i,2]))
  }
  
  IMP <- result$Impute_All
  IMP[posirows,]<-0
  
  PTSZ=list()
  PTDO=list()
  for(j in 1:ncol(observed)){
    indexobserved0=(observed[,j]==0)
    predictv=IMP[,j][indexobserved0]
    Truevalue=expected[,j][indexobserved0]
    
    PTSZ[[j]]=sum(Truevalue==0 & predictv==0)/sum(Truevalue==0)
    PTDO[[j]]=sum(Truevalue>0 & predictv>0)/sum(Truevalue>0)
  }
  
  PTSZ=unlist(PTSZ)
  PTDO=unlist(PTDO)
  summa_mean=data.frame(PTSZ=mean(PTSZ), SD1=sd(PTSZ), thresh=tt)
  rownames(summa_mean)=NULL
  return(summa_mean)
}
