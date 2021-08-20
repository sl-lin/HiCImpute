#' 
#' This function analyzes both simulated and real datasets, depending on the inputs of the functions.
#'
#' @param result Output of MCMCImpute for simulated data or the organized results of real data.
#' @param cell_index Indicates which cell is used to draw heatmaps and scatterplot.
#' @param n Dimension of 2D contact matrix.
#' @param cell_type A vector of underlying true cluster.
#' @param dims The dimension of 2D matrix.
#' @param perplexity numeric; Perplexity parameter (should not be bigger than \eqn{3\times perplexity<nrow(X)-1}).
#' @param seed Random seed for generating t-SNE data.
#' @param kmeans Logical, whether apply K-means clustering on the t-SNE data.
#' @param ncenters Number of centers in K-means clustering analysis.
#'
#' @return A list of accuracy measurements and plots.
#' @export
#'
#' @import gridExtra
#'
#' @examples
#' data("K562_1_true")
#' options(digits = 2)
#' scHiC_assess(result=K562_T1_4k_result)

scHiC_assess <- function(result, cell_index=1, n, cell_type, dims = 2,perplexity=10, seed=1000, 
                         kmeans = TRUE, ncenters = 2){

   observed = result$scHiC
   expected = result$Expected
  
  if(is.null(expected)){##### for real data analysis
    
    imputed=result$Imputed
    
    ############# correlation
    corr=NULL
    for (i in 1:ncol(observed)) {
      ind=observed[,i]!=0
      corr[i]=cor(observed[,i][ind],imputed[,i][ind])
    }
    corr=as.data.frame(corr)
    p1=ggplot(corr, aes(x="corr", y=corr))+geom_boxplot()+ggtitle("Boxplot of correlations") + theme(axis.title.x=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank()) 
    
    ############# scatterplot
    p2=SOVI(obsvec=observed[,cell_index], impvec = imputed[,cell_index])
    
    p3=scHiC_tSNE(observed, cell_type=cell_type, 
                  dims = dims, perplexity=perplexity, seed=seed, title="Observed", 
                  kmeans = kmeans, ncenters = ncenters)
    
    p4=scHiC_tSNE(imputed, cell_type=cell_type, 
                  dims = dims, perplexity=perplexity, seed=seed, title="HiCImpute", 
                  kmeans = kmeans, ncenters = ncenters)
     
    p=ggarrange(p1, p2, p3, p4,
              labels = c("A", "B", "C", "D"),
              ncol = 2, nrow = 2)
    
    
    cluster1 = scHiC_Kmeans(observed, centers=2, nstart=100, iter.max=1000, seed=1)
    dd1 = as.data.frame(cbind(cluster=cluster1$cluster,type=cell_type))
    table1 = table(dd1$cluster,dd1$type)
    
    cluster2 = scHiC_Kmeans(imputed, centers=2, nstart=100, iter.max=1000, seed=1)
    dd2 = as.data.frame(cbind(cluster=cluster2$cluster,type=cell_type))
    table2 = table(dd2$cluster,dd2$type)
    
    mylist=list(p, table1, table2)
    names(mylist)=c("plots", "obs_cluster", "imp_cluster")
    return(mylist)
  }else{ ######### for simulated data
    
    ##readin data of different types
    
    if(is.list(expected)){
      EXP=NULL
      for (c in 1:length(expected)) {
        EXP=cbind(EXP, mattovec(expected[[c]]))
      }
    }else if(is.array(expected) & length(dim(expected))>2){
      EXP=NULL
      for (c in 1:(dim(expected)[3])) {
        EXP=cbind(EXP, mattovec(expected[,,c]))
      }
    }else if(is.matrix(expected)){
      EXP=expected
    }else{
      message("The input type must be a matrix, 3D array or a list of 2D matrix.")
      exit()
    }
    
    expected = EXP
    
    imputed = result$Impute_SZ
      
    PTSZ=list()
    PTDO=list()
    #MSE=list()
    CIEZ=list()
    CIEA=list()
    AEOZ=data.frame()
    AEOA=data.frame()
    
    for(j in 1:ncol(observed)){
      index0 <- observed[,j]==0
      index1 <- observed[,j]==1
      index2 <- observed[,j]>=2
      
      ## calibration
      m1 <- lm(result$Impute_SZ[,j][observed[,j]>0]~observed[,j][observed[,j]>0])
      impu0 <- result$Impute_SZ[,j][result$Impute_SZ[,j]>0 & observed[,j]==0]
      obs2 <- (impu0-m1$coefficients[1])/(m1$coefficients[2])
      result$Impute_SZ[,j][result$Impute_SZ[,j]>0 & observed[,j]==0] <- obs2
      
      indexnotrue0=(expected[,j]!=0)
      
      predictv=result$Impute_SZ[,j][index0]
      Truevalue=expected[,j][index0]
      
      PTSZ[[j]]=sum(Truevalue==0 & predictv==0)/sum(Truevalue==0)
      PTDO[[j]]=sum(Truevalue>0 & predictv>0)/sum(Truevalue>0)
      
      #MSE[[j]]=mean((result$Impute_SZ[,j]-expected[,j])^2)
      
      AEOZ=rbind(AEOZ, data.frame(AE=(abs(predictv-Truevalue)),Sample=j))
      AEOA=rbind(AEOA, data.frame(AE=abs(result$Impute_SZ[,j]-expected[,j]),Sample=j))
      
      CIEZ[[j]]=cor(Truevalue,predictv)
      CIEA[[j]]=cor(result$Impute_SZ[,j],observed[,j])
    }
    
    PTSZ=unlist(PTSZ)
    PTDO=unlist(PTDO)
    #MSE=unlist(MSE)
    CIEZ=unlist(CIEZ)
    CIEA=unlist(CIEA)
    
    summa_mean=data.frame(PTSZ=mean(PTSZ),PTDO=mean(PTDO),
                          AEOZ=mean(AEOZ$AE),AEOA=mean(AEOA$AE),
                          CIEZ=mean(CIEZ),CIEA=mean(CIEA))
    summa_se=data.frame(PTSZ=sd(PTSZ),PTDO=sd(PTDO),
                        AEOZ=sd(AEOZ$AE),AEOA=sd(AEOA$AE),
                        CIEZ=sd(CIEZ),CIEA=sd(CIEA))
    rownames(summa_mean)=NULL
    rownames(summa_se)=NULL
    
    
    ############# fix PTSZ=0.95
    PTSZ_95 = PTSZ95(observed=observed, expected=expected, result=result)
    
    par(mar = c(3,3,3,5))
    par(mfrow=c(2,3))
    scHiC_hm(observed[,cell_index], n, title="Observed")
    fig_label("A", cex=2) 
    scHiC_hm(expected[,cell_index], n, title="Expected")
    fig_label("B", cex=2) 
    scHiC_hm(result$Impute_SZ[,cell_index], n, title="HiCImpute")
    fig_label("C", cex=2) 
    ############# ROC curve
    scHiC_ROC(observed=observed, expected=expected, result=result)
    fig_label("D", cex=2) 
    SEVI(obsvec=observed[,cell_index], expvec=expected[,cell_index], impvec=result$Impute_All[,cell_index] )
    fig_label("E", cex=2) 

    mylist <- list(summa_mean, summa_se, PTSZ_95)
    names(mylist) <- c("summary_mean", "summary_se", "PTSZ_95")
    return(mylist)
  }
}


