#' 
#' This is the flagship function of HiCImpute. It identifies structural zeros and imputes 
#' sampling zeros under a Bayesian framework. The outputs can be used to facilitate downstream 
#' analysis such as clustering, subtype discovery, or 3D structure construction.
#' 
#' @param scHiC The single-cell Hi-C matrix. It can take three types of formats. The preferred format is a single-cell matrix with 
#' each column being a vector of the upper triangular matrix without including the diagonal entries 
#' of the 2D matrix of a single-cell. Another types of formats are a list with each element being a 
#' 2D single-cell contact matrix, or a 3D (\eqn{n\times n\times k}) array that has k matrices of dimension 
#' \eqn{n\times n}. HiCImpute automatically transforms these two types of input into a matrix with each 
#' column being the vector of upper triangular matrix of a single-cell. For a single-cell matrix of 
#' size \eqn{n \times n}, the length of the vector should be \eqn{n\times(n-1)/2}. We only need the upper 
#' triangular matrix because the Hi-C matrix are symmetrical. 

#' @param bulk The bulk data. It can take two types of formats. A 2D bulk matrix of dimension \eqn{n\times n} 
#' or a vector of the upper triangular entries of 2D bulk matrix. It can provide information for priors settings. 
#' If bulk data is not available, simply set it to be NULL, and MCMCImpute will sum up the single-cells to 
#' construct a bulk data.
#'
#' @param expected Underline true counts of the simulated data. For real data analysis, just set it as NULL.It takes
#'three formats that is the same as scHiC.
#' @param startval The starting value for the vector of parameters 
#' \eqn{\Theta=(\alpha, \mu^\gamma, \beta, \mu, a, \delta, b, \pi, s, \mu_1, \cdots, \mu_K)}. 
#' See xie et al. for the details of these parameters. The default value is as set in the function.
#' @param n Integer. The dimension of single-cell matrix.
#' @param epsilon1 The range size of  \eqn{\delta} that is used to monitor the prior mean of \eqn{\pi_{ij}}, 
#' the probability that the pair \eqn{(i,j)} do not interact. The default value of \eqn{\epsilon_1} is 0.5.
#' @param epsilon2 The range size of \eqn{B} that is used to monitor the prior mean of \eqn{\mu_{ij}}, the intensity 
#' of interaction between pair \eqn{(i,j)}. The default value  of \eqn{\epsilon_2} is 5.
#' @param mc.cores The number of cores to be used in mclapply function that can parallelly impute the 
#' matrix. The default value is 1 (no parallelization), but the users is advised to a higher number to increase 
#' computational speed if their computer has parallel computing capability.
#' @param cutoff The threshold of \eqn{\pi_{ij}} that is used to define structural zeros. The default value is 0.5. 
#' That is, if the probability of being a SZ is greater than 0.5, then the pair \eqn{(i,j)} are labelled as not 
#' interacting due to underlying biological mechanism.
#' @param niter  The number of iterations for the MCMC run. Default is 30000.
#' @param burnin The number of burn-in iteration. Default is 5000.
#' 
#' @import parallel
#' @import Rcpp
#' @import RcppArmadillo
#' 
#' @return A list of posterior mean of probability (SZ), the imputed data without defining SZ (Impute_All), 
#' and imputed data with SZ, using the threshold (Impute_SZ). 
#' @export
#'
#' @examples
#' data("K562_T1_4k")
#' data("K562_bulk")
#' data("K562_T1_4k_true")
#' scHiC=K562_T1_7k
#' set.seed(1234)
#' T1_4k_result1=MCMCImpute(scHiC=K562_T1_4k,bulk=K562_bulk, expected=K562_T1_4k_true,
#'                          startval=c(100,100,10,8,10,0.1,900,0.2,0,replicate(dim(scHiC)[2],8)),n=61,
#'                                                   mc.cores = 1,cutoff=0.5, niter=50,burnin=10)

MCMCImpute <- function(scHiC, bulk=bulk, expected=NULL,
                       startval=c(100,100,10,8,10,0.1,900,0.2,0,
          replicate(dim(single)[2],8)),n,epsilon1=0.5,epsilon2=5,mc.cores = 1,cutoff=0.5, niter=30000,burnin=5000){
  ##readin data of different types
  if(is.list(scHiC)){
    single=NULL
    for (c in 1:length(scHiC)) {
      single=cbind(single, mattovec(scHiC[[c]]))
    }
  }else if(is.array(scHiC) & length(dim(scHiC))>2){
    single=NULL
    for (c in 1:(dim(scHiC)[3])) {
      single=cbind(single, mattovec(scHiC[,,c]))
    }
  }else if(is.matrix(scHiC)){
    single=scHiC
  }else{
    stop("The input type must be a matrix, 3D array or a list of 2D matrix.")
  }
  
  if (!is.vector(bulk)){
    bulk = mattovec(bulk)
  }
    
  if (nrow(single) != length(bulk)){
    stop("The length of single cells and bulk data are different!")
  }
  
  if (nrow(single) != (n*(n-1)/2)){
    stop("The provided n and single cells do not match in size!")
  }
  
  ## normalize single cells to get same depth as max
  single_sum <- apply(single, 2, sum)
  max_single <- max(single_sum)
  lambda <- single_sum/max_single
  single_norm <- t(t(single)/lambda)
  
  ## variance for each contact pair
  mu_sigma = neivar(single_norm, nei=5, n)
  mu_sigma_vec = mu_sigma[upper.tri(mu_sigma)]
  
  ## single weighted sum 
  n_single = dim(single)[2]
  weight <- single_sum/sum(single_sum)
  type_bulk_vec <- apply(t(t(single_norm)*weight),1,sum)
  type_bulk <- matrix(0, n, n) 
  type_bulk[upper.tri(type_bulk, diag=FALSE)] <- type_bulk_vec
  
  ## correction matrix
  correct = correctfac(type_bulk,nei = 5)
  if(is.null(bulk)){
    bulk =  apply(single, 1, sum)
  }
  
  ## combine information into one matrix
  m <- single
  rowmax <- apply(m, 1, max) 
  prop0 <- rowMeans(m==0, na.rm = FALSE, dims = 1)
  correctvec <- correct[upper.tri(correct,diag = FALSE)]
  B <- bulk*correctvec*rowmax/sum(bulk) 
  m <- as.matrix(cbind(m,rowmax,prop0,correctvec,B))
  
  ## apply oneimpute function to each row of matrix
  result=mclapply(0:(nrow(single)-1), oneimpute, niter=niter, burnin=burnin, n_single=n_single, m=as.matrix(m), lambda=lambda, mu_sigma_vec=mu_sigma_vec, startval=startval, n=n, epsilon1=epsilon1, epsilon2=epsilon2, mc.cores = mc.cores)
  result=do.call(rbind,result)
  
  ## organize MCMC result
  pii_mean <- result[,1]
  pii <- matrix(0, n, n)
  pii[upper.tri(pii, diag=FALSE)] <- pii_mean
  post_mu_k <- result[,-1]
  
  ## true 0 positions figured out by MCMC
  # posi <- which(pii>=cutoff,TRUE)
  # posirows <- NULL
  # for (i in 1:(dim(posi)[1])) {
  #   posirows <- c(posirows,matrow(posi[i,1],posi[i,2]))
  # }
  
  ## Imputation
  IMP1 <- NULL
  for (i in 1:ncol(single)) {
    imp <- post_mu_k[,i]*lambda[i]
    IMP1 <- cbind(IMP1,imp)
  }
  
  IMP2=IMP1
  IMP2[pii_mean>0.5,]<-0
  
  ## output
  output=list(single,pii,IMP1, IMP2, expected)
  names(output) = c("scHiC","SZ","Impute_All","Impute_SZ", "Expected")
  return(output)
}



