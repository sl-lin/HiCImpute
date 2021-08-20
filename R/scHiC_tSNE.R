#' This function visualize scHi-C data using t-SNE (t-distributed stochastic neighbor embedding) 
#' and applying Kmeans clustering followed by xie et al. 2021.
#'
#' @param data The observed matrix, with each column being the uppertriangular of a single cell HiC matrix.
#' @param cell_type A vector that indicates cell type.
#' @param dims Integer. Output dimentionality. Default=2.
#' @param perplexity Numeric; Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation).
#' @param check_duplicates Logical; Checks whether duplicates are present. It is best to make sure there are no duplicates present and set this option to FALSE, especially for large datasets (default: TRUE).
#' @param seed Random seed.
#' @param ncenters Number of clusters in kmeans clustering.
#' @param title Title of the plot.
#'
#' @return
#' A stne visualization plot.
#' @export
#'
#' @import Rtsne
#'
#' @import ggpubr
#' 
#' @import ggplot2
#'
#' @examples
#' scHiC_tSNE(GSE117874_chr1_wo_diag, cell_type=c(rep("GM",14),rep("PBMC",18)), 
#' dims = 2,perplexity=10, seed=1000, title="Observed GSE117874", 
#' kmeans = TRUE, ncenters = 2)

scHiC_tSNE <- function(data, cell_type, dims = 2, perplexity=10, check_duplicates = FALSE, seed=1234, title=NULL,
                      kmeans=TRUE, ncenters) {
  
  mydata <- t(scale(data)) # standardize variables
  #rownames(mydata) = c(paste("A",1:131,sep="_"), paste("B",1:180,sep="_"))
  
  set.seed(seed)
  tsne_dat <- Rtsne(mydata, perplexity=perplexity, check_duplicates = check_duplicates)
  
  data=tsne_dat
  if(kmeans){
    my.xy = data.frame(x=data$Y[,1], y=data$Y[,2])
    
    res.km <- kmeans(my.xy, centers = ncenters, nstart = 1)
    
    my.xy$cluster <- factor(res.km$cluster)
    my.xy$type <- cell_type
    
    ggscatter(
      my.xy, x = "x", y = "y", xlab = "",
      ylab = "", title=title, color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
      shape = "type", size = 2,  legend = "right", ggtheme = theme_bw()) +
      scale_shape_manual(values=c(2, 19, 17, 10)[1:ncenters])
  }else{
    my.xy = data.frame(x=data$Y[,1], y=data$Y[,2], group=cell_type)
    ggplot(my.xy) + geom_point(aes(x=x, y=y, color=group))
  }
}






