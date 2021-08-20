#' This function generates scatterplot of expected versus imputed. 
#'
#' @param obsvec A vector of observed single cell.
#' @param expvec A vector of expected single cell.
#' @param impvec A vector of imputed single cell.
#'
#' @return The scatterplot of expected versus imputed, with read dots being the observed zero pairs.
#' @export
#'
#' @examples
#' SEVI(obsvec=K562_T1_7k[,1], expvec=K562_1_true[,1], impvec=T1_7k_imp[,1] )

SEVI <- function(obsvec, expvec, impvec) {
   plot(expvec, impvec, pch=20, cex=0.4, col="blue", xlab = "Expected", ylab = "Imputed",main = "SEVI")
   points(expvec[obsvec==0], impvec[obsvec==0], pch=20, cex=0.4, col="red")
}

