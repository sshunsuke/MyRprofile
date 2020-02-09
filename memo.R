
# draw lines/points. 
#   plotMats(mat1, mat2, ...)    # (length(mat1) == length(mat2)) == FALSE
# If data have same length of rows, you should use matplot(). 
plotMats <- function(..., xlab='X', ylab='Y', main=NULL, type='l', log='',
                     xlim=NULL, ylim=NULL, x0=TRUE, y0=TRUE) {
  mats <- list(...)
  plotMats2(mats, xlab=xlab, ylab=ylab, main=main, type=type, log=log,
            xlim=xlim, ylim=ylim, x0=x0, y0=y0)
}

plotMats2 <- function(mats, xlab='X', ylab='Y', main=NULL, type='l', log='',
                      xlim=NULL, ylim=NULL, x0=TRUE, y0=TRUE) {
  len <- length(mats)
  df <- switch(type, 'p' = points, 'l' = lines)
  
  if (is.null(xlim) && length(grep('.*x.*', log)) > 0) {
    if (x0) { x0 <- FALSE; cat("x0 is forcedly turned off for log axis.") }
  }
  if (is.null(ylim) && length( grep('.*y.*', log) ) > 0) {
    if (y0) { y0 <- FALSE; cat("y0 is forcedly turned off for log axis.") }
  }
  
  rangeX <- ifelse(x0, 0, mean( (mats[[1]])[,1] ))
  rangeY <- ifelse(y0, 0, mean( (mats[[1]])[,1] ))
  for (i in 1:len) {
    rangeX <- range( rangeX, (mats[[i]])[,1] )
    rangeY <- range( rangeY, (mats[[i]])[,2] )
  }
  if (is.null(xlim)) { xlim <- rangeX }
  if (is.null(ylim)) { ylim <- rangeY }
  
  plot(mats[[1]], xlab=xlab, ylab=ylab, main=main, type=type, xlim=xlim, ylim=ylim, log=log)
  if (len > 1) {
    for (i in 2:len) {
      df(mats[[i]], col=i)
    }
  }
}






testNewtonMethod <- function(roughness, D, Re) {
  fun <- function(fD) {
    (1 / sqrt(fD)) + 2 * log10( roughness / D / 3.71 + 2.51 / Re / sqrt(fD))
  }
  
  # Derivative of fun().
  dFun <- function(fD) {
    - fD^(-3/2) * (1/2 + 2.51 / log(10) / (2.51 / Re / sqrt(fD) + roughness / 3.71 / D) / Re)
  }
  
  fD_0 <- FCP$fD.Blasius(Re)
  newtonRaphson(fun, dFun, fD_0)
}


testBisection <- function(roughness, D, Re) {
  fun <- function(fD) {
    (1 / sqrt(fD)) + 2 * log10( roughness / D / 3.71 + 2.51 / Re / sqrt(fD))
  }
  
  fD_0 <- FCP$fD.Blasius(Re)
  bisection(fun, fD_0/2, fD_0*2)
}


FCP$fD.Colebrook(0.0005, 0.02, 8000)
testNewtonMethod(0.0005, 0.02, 8000)
testBisection(0.0005, 0.02, 8000)
