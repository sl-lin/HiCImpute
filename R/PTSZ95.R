#'This function calculates PTDO when fix PTSZ=0.95.
#' @param observed Observed single cells matrix with each column being the upper triangular of a single cell.
#' @param expected Underline true counts from simulation.
#' @param result Result form MCMCImpute function.
#'
#' @return A vector of PTDO and its SD when fixing PTSZ to be 0.95, and the threshold used in that case.
#' @export
#'
#' @examples
#' PTSZ95(observed=K562_T1_7k, expected=K562_1_true, result=T1_7k_res)

PTSZ95 <- function(observed, expected, result){
  
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
    return(mean(PTSZ))
  }
  
  threshold <- seq(0.001, max(result$SZ), length.out = 200)
  PTSZout <- unlist(lapply(threshold,PTSZfun))
  tt <- rev(threshold)[which.min(abs(rev(PTSZout)-0.95))]
  
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
  summa_mean=data.frame(PTDO=mean(PTDO), SD2=sd(PTDO), thresh=tt)
  rownames(summa_mean)=NULL
  return(summa_mean)
}
