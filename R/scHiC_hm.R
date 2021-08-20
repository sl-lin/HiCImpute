#' 
#' This function draws heatmap of HiC data so that we can visually compares the imputation results.
#'
#' @param datvec A vector of upper triangular mamtrix.
#' @param n Dimension of 2D matrix (i.e., the number of segments).
#' @param title The title of the heatmap.
#'
#' @return Heatmap of the matrix.
#' @export
#'
#' @examples.
#' data("K562_1_true")
#' scHiC_hm(K562_1_true[,1], 61, title="Expected")
scHiC_hm <- function(datvec, n, title="Heatmap") {
  library(plsgenomics)
  normmatrix <- function(matr) {
    maxvalue <- max(matr[upper.tri(matr)])
    minvalue <- min(matr[upper.tri(matr)])
    normmatr <- (matr-minvalue)/(maxvalue-minvalue) 
    return(normmatr)
  }
  
  mat <- matrix(0, n, n)
  mat[upper.tri(mat, diag=FALSE)] <- datvec
  return(matrix.heatmap(normmatrix(mat), main=title))
}

