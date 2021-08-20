#'This function simulates single cells from 3D structure.
#' @param data 3D coordinates of single cell.
#' @param alpha_0 Parameter that controls sequence depth of data.
#' @param alpha_1 Parameter that controls sequence depth of data.
#' @param beta_l Parameter that controls effect size of covariate.
#' @param beta_g Parameter that controls effect size of covariate.
#' @param beta_m Parameter that controls effect size of covariate.
#' @param gamma Quantile that is used as the threshold.
#' @param eta Percent of structural zeros that are set to be common structural zeros among all single-cells.
#' @param n_single Number of single cells to be generated.
#'
#' @return A list of underline true count, SZ positions, and generated single cells.
#' @export
#'
#' @examples
#' #Load 3d structure generated from SIMBA package
#' load("simba_3strs.rdata")
#' Set random seed
#' set.seed(1234)
#' #Generate 100 random type1 single cells
#' simudat <- scHiC_simulate(data=str1, alpha_0=5.6,alpha_1=-1, beta_l=0.9,beta_g=0.9,
#' beta_m=0.9,gamma=0.1,eta=0.8, n_single=10) 

scHiC_simulate <- function(data=str1,alpha_0,alpha_1, beta_l,beta_g,beta_m,gamma,eta,n_single){  
  
##define function that summarize information from 3d coordinates.
    nposi <- dim(data)[1]
    
    position <- cbind(data,runif(nposi,min = 0.2, max = 0.3),runif(nposi,min = 0.4, max = 0.5),runif(nposi,min = 0.9, max = 1))
    
    ##matrix containing distance info
    distance <- matrix(0,nrow=nposi,ncol=nposi)
    for (i in 1:(nposi-1)) {
      for (j in (i+1):nposi) {
        distance[i,j] = sqrt(sum((position[i,1:3]-position[j,1:3])^2)) 
      }
    }
    
    ##matrix containing lambda info
    lambda <- matrix(0,nrow=nposi,ncol=nposi)
    for (i in 1:(nposi-1)) {
      for (j in (i+1):nposi) {
        l = position[i,4]*position[j,4]
        g = position[i,5]*position[j,5]
        m = position[i,6]*position[j,6]
        lambda[i,j] <- exp(alpha_0+alpha_1*log(distance[i,j])+beta_l*log(l)+beta_g*log(g)+beta_m*log(m))
      }
    }
    
    ##sequence depth
    seqdepth <- sum(lambda[upper.tri(lambda)])
    ##true 0 positions
    thresh <- quantile(lambda[upper.tri(lambda)],eta)
    posi0 <- which(lambda<thresh & upper.tri(lambda),TRUE)
    true0 <- posi0[sample(nrow(posi0),size=0.5*nrow(posi0),replace=FALSE),]  
  
  ## Define function
  rpoi <- function(x) {  ## function to generate one Poisson random number
    r <- rpois(1,x)
    return(r)
  }
  
  matrow <- function(x,y) {   ## function to calculate rows corresponding to location
    r <- x + (y-1)*(y-2)/2
    return(r)
  }
  
  # function to downsample a matrix
  downsamplemat <- function(mat, samplerate=0.5) {
    new <- matrix(0, nrow(mat), ncol(mat))
    #colnames(new) <- colnames(mat)
    #rownames(new) <- rownames(mat)
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        new[i,j] <- sum(runif(mat[i,j], 0, 1) < samplerate)
      }
    }
    return(new)
  }
  
  # v is the vector of full sample
  subsampling=function(v,eta)
  {
    # set.seed(1234)
    out=rep(0,length(v))
    v1=v
    # 0 should also include in the sampling
    v1[v==0]=1
    v2=cumsum(v1)
    # start
    v3=cumsum(v1)-v+1
    total=sum(v1)
    index=sample(1:total,total*eta,replace=F)
    for(i in 1:length(index))
    {
      # ith number
      intersect.row=(v3<=index[i] & v2>=index[i])
      out[intersect.row]=out[intersect.row]+1
    }
    # where input is 0, output should also be 0
    out[v==0]=0
    return(out)
  }
  # check
  #v=rpois(1000,2)
  #test=subsampling(v,0.2)
  #sum(test)/sum(v)
  #sum(test==0)/sum(v==0)
  #sum(test<=v)
  
  l <- dim(true0)[1]
  
  #random 0 positions
  if(eta==0){
    #underline true model for all single cells  
    underline <-  lambda
    underline <- underline[upper.tri(underline,diag = FALSE)]
    single <- matrix(rep(underline,each=n_single), ncol=n_single, byrow=TRUE)
    random0 <- NULL
    my_list <- list(single, random0)
    names(my_list) <- c("single","random0")
    return(my_list)
  }else{
    true0rows <- NULL
    for (i in 1:(dim( true0)[1])) {
      true0rows <- c(true0rows,matrow( true0[i,1], true0[i,2]))
    }
    
    random0 <- true0rows[sample(1:l,floor(l*eta))]
    
    #underline true model for all single cells  
    underline <-  lambda
    underline <- underline[upper.tri(underline,diag = FALSE)]
    randlam <- underline[random0]
    underline[true0rows] <- 0
    
    truecount <- NULL
    for (j in 1:n_single) {
      s <- underline
      s[random0] <- randlam*rbinom(length(random0),1,0.5)
      truecount <- cbind(truecount,s)
    }
    
    singledat = apply(truecount,c(1,2),rpoi)
    
    my_list <- list(truecount, random0, singledat)
    names(my_list) <- c("truecount","random0", "singledat")
    return(my_list)
  }
}


