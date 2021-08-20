correctfac <- function(type_bulk,nei){
  n = dim(type_bulk)[1]
  mean_bulk = sum(type_bulk[upper.tri(type_bulk)])/sum(upper.tri(type_bulk))
  correct = matrix(rep(0,n*n),nrow=n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      m <- type_bulk[max(1,i-nei):min(n,i+nei), max(1,j-nei):min(n,j+nei)]
      c <- mean(m[!is.na(m)])
      correct[i,j] = c
    }
  }
  return(correct)
}


matrow <- function(x,y) {   
  r <- x + (y-1)*(y-2)/2
  return(r)
}


mattovec <- function(mat){
  vec <- mat[upper.tri(mat, diag = FALSE)]
  return(vec)
}


neivar <- function(single, nei, n){
  single_mat <- list()
  for (k in 1:dim(single)[2]) {
    m <- matrix(NA,n,n)
    m[upper.tri(m)] <- single[,k]
    single_mat[[k]] <- m
  }
  correct = matrix(0 ,nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      neighbor <- NULL
      for (l in 1:dim(single)[2]) {
        neighbor = c(neighbor,as.vector(single_mat[[l]][max(1,i-nei):min(n,i+nei), max(1,j-nei):min(n,j+nei)]))
      }
      c <- neighbor[!is.na(neighbor)]
      cc <- c[c>0]
      correct[i,j] = sd(cc)
    }
  }
  return(correct)
}




fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

