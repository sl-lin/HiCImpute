#'This package draws ROC (Receiver operating characteristic) curve to visually demonstrate ability 
#'to tell SZ from DO. 
#'
#' @param observed Observed sngle cell with each column being the upper triangular of single cell.
#' @param expected Underline true count of simulated data.
#' @param result Result from MCMCImpute funciton.
#'
#' @return A plot of ROC curve.
#' @export
#'
#' @import DescTools
#' 
#' @examples
#' scHiC_ROC(observed=K562_T1_7k, expected=K562_1_true, result=T1_7k_res)

scHiC_ROC <- function(observed, expected, result){
  
  observed_sum <- apply(observed, 2, sum)
  max_observed <- max(observed_sum)
  observedlam <- observed_sum/max_observed
  
  ## true 0 positions figured out by MCMC
  PTSZfun <- function(thresh) {
    posi <- which(result$SZ>=thresh,TRUE)
    posirows <- NULL
    for (i in 1:(dim(posi)[1])) {
      posirows <- c(posirows,matrow(posi[i,1],posi[i,2]))
    }
    IMP <-result$Impute_All
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
    return(mean(PTSZ))
  }
  
  threshold <- seq(0.001, max(result$SZ), length.out = 200)
  PTSZout <- unlist(lapply(threshold,PTSZfun))
  
  ptsz_range <- seq(0,1,length.out = 100)
  output=NULL
  for (k in 1:length(ptsz_range)) {
    tt <- rev(threshold)[which.min(abs(rev(PTSZout)-ptsz_range[k]))]
    
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
    summa_mean=data.frame(PTSZ=mean(PTSZ),PTDO=mean(PTDO),thresh=tt)
    output <- rbind(output,summa_mean)
  }
  auc=format(round(AUC(c(1-output$PTSZ,0),c(output$PTDO,0)),2), nsmall = 2)
  
  #return(output)
  plot(c(1-output$PTSZ,0),c(output$PTDO,0), type = "l", xlab = "1-PTSZ", ylab = "PTDO",
      main = paste("ROC, AUC=",auc,sep = ""),  pch=3,col="darkorchid", lwd=2, ylim = c(0,1))
}

