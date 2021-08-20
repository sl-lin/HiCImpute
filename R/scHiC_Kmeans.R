#' This function conduct Kmeans clustering analysis on scHi-C data. 
#'
#' @param data The observed or imputed matirx, with each column being the uppertriangular of a single cell HiC matrix.
#' @param centers Either the number of clusters, say k, or a set of initial (distinct) cluster centres. If a number, a random set of (distinct) rows in x is chosen as the initial centres.
#' @param nstart If centers is a number, how many random sets should be chosen.
#' @param iter.max The maximum number of iterations allowed.
#' @param seed Random seed.
#' @param algorithm Character: may be abbreviated. Note that "Lloyd" and "Forgy" are alternative names for one algorithm.
#' @param trace Logical or integer number, currently only used in the default method ("Hartigan-Wong"): if positive (or true), tracing information on the progress of the algorithm is produced. Higher values may produce more tracing information.
#'
#' @return Kmeans clustering results.
#' @export
#'
#' @examples
#' data("GSE117874_chr1_wo_diag")
#' data("GSE117874_imp")
#' cluster=scHiC_Kmeans(GSE117874_chr1_wo_diag, centers=2, nstart=1, iter.max=1000, seed=1)

scHiC_Kmeans <- function(data, centers, nstart=50, iter.max=200, seed=1234,
                           algorithm = c("Hartigan-Wong", "Lloyd", "Forgy",
                                         "MacQueen"), trace=FALSE) {
  mydata <- t(scale(data)) # standardize variables
  set.seed(seed)
  CLUSTER <- kmeans(mydata, centers = centers, nstart=nstart ,iter.max=iter.max,
                    algorithm = algorithm, trace=FALSE)
  #adjustedRandIndex(CLUSTER$cluster, c(rep(1,131),rep(2,180)))
  return(CLUSTER)
}




