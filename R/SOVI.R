#' This function generates scatterplot of observed versus imputed for nonzero observed counts.
#'
#' @param obsvec A vector of observed single cell.
#' @param impvec A vector of imputed single cell.
#'
#' @return The scatterplot of observed versus imputed.
#' @export
#' 
#' @import ggplot2
#'
#' @examples 
#' data("GSE117874_imp")
#' data("GSE117874_chr1_wo_diag")
#' SOVI(obsvec = GSE117874_chr1_wo_diag[,1], impvec = GSE117874_imp[,1])

SOVI <- function(obsvec, impvec) {
  dat=as.data.frame(cbind(obs=obsvec[obsvec>0], imp=impvec[obsvec>0]))
  ggplot(dat, aes(x=obs, y=imp)) + geom_point() + ggtitle("SOVI") + xlab("Observed") + ylab("Imputed")
  #plot(observed[observed>0], imputed[observed>0], pch=20, cex=0.4, col="blue", xlab = "Observed", ylab = "Imputed",main = "SOVI")
}