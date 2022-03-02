
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





# -----------------------------------------------------------------------------
# Properties of Methane
# 
# Reference:
#   Duan, Z., Moller, N., & Weare, J. (1992). An Equation of State for the 
#   CH4-CO2-H2O System: I. Pure systems from 0-Degrees-C to 1000-Degrees-C and 
#   0 to 8000 Bar. Geochimica Et Cosmochimica Acta, 56(7).
#   doi:http://dx.doi.org/10.1016/0016-7037(92)90347-L
#
#   Friend, D.G., Ely, J.F., Ingham, H. (1989). Thermophysical Properties of
#   Methane. J. Phys. Chem. Ref. Data 18, 583.
# -----------------------------------------------------------------------------

# Z Factor of Methane (Compressibility Factor), Duan (1992)
methaneZ <- function(P, T, range=c(1e-8,1000), tol=1e-08) {
  a1  <-  8.72553928e-02
  a2  <- -7.52599476e-01
  a3  <-  3.75419887e-01
  a4  <-  1.07291342e-02
  a5  <-  5.49626360e-03
  a6  <- -1.84772802e-02
  a7  <-  3.18993183e-04
  a8  <-  2.11079375e-04
  a9  <-  2.01682801e-05
  a10 <- -1.65606189e-05
  a11 <-  1.19614546e-04
  a12 <- -1.08087289e-04
  
  alpha <- 4.48262295e-02
  beta  <- 7.5397e-01
  gamma <- 7.7167e-02
  
  Pc <- 46.41 * 100 * 1000   #  46.71 (bar)
  Tc <- -82.55 + 273.15      # -82.55 (degC)
  Vc <- R * Tc / Pc
  
  Pr = P / Pc   # Eq (3)
  Tr = T / Tc   # Eq (4)
  
  # Calculate 'Vr' by Newton-Raphson method - Eq (1)
  fVr <- function(Vr) {
    B <- a1 + (a2 / Tr^2) + (a3 / Tr^3)     # Eq (2a)
    C <- a4 + (a5 / Tr^2) + (a6 / Tr^3)     # Eq (2b)
    D <- a7 + (a8 / Tr^2) + (a9 / Tr^3)     # Eq (2c)
    E <- a10 + (a11 / Tr^2) + (a12 / Tr^3)  # Eq (2d)
    F <- alpha / Tr^3                       # Eq (2e)
    
    l  <- (Pr / Tr) * Vr
    r1 <- 1 + B/Vr + C/Vr^2 + D/Vr^4 + E/Vr^5
    r2 <- (F/Vr^2) * (beta + gamma / Vr^2) * exp(-gamma / Vr^2)
    
    l - r1 - r2
  }
  f <- uniroot(fVr, range, tol=tol)
  Vr = f$root
  
  # DEBUG( "P:  %.3f(kPa), T: %.3f(K)\n", P/1000, T )
  # DEBUG( "Pr: %f, Tr: %f, Vr: %f\n", Pr, Tr, Vr )
  # print(f)
  
  Pr * Vr / Tr
}

# Viscosity of methane, Friend & Ingham (1989)
methaneViscosity <- (function(){
  vis0 <- function(T) {
    c_ <- c( -3.0328138281,
             16.918880086,
             -37.189364917,
             41.288861858,
             -24.615921140,
             8.9488430959,
             -1.8739245042,
             0.20966101390,
             -9.6570437074*(10^-3) )
    t_ <- T / 174
    f <- function(i) { c_[i] * (t_ ^ ((i-1)/3 - 1)) }
    omega <- ( sum(mapply(f, 1:9)) )^(-1)
    10.5 * (t_^0.5) / omega
  }
  
  visEx <- function(density, T) {
    r_ <- c(1, 1, 2, 2,   2, 3, 3, 4, 4, 1, 1)
    s_ <- c(0, 1, 0, 1, 1.5, 0, 2, 0, 1, 0, 1)
    g_ <- c( 0.41250137,
             -0.14390912,
             0.10366993,
             0.40287464,
             -0.24903524,
             -0.12953131,
             0.06575776,
             0.02566628,
             -0.03716526,
             -0.38798341,
             0.03533815)
    
    rDensity <- (density/methaneM) / (10.139 * 10^3)  # reduced density
    rT <- T / 190.551             # reduced temeperature
    f <- function(i) { g_[i] * rDensity^(r_[i]) * rT^(s_[i]) }
    12.149 * sum(mapply(f, 1:9)) * (1 + sum(mapply(f, 10:11)))^(-1)
  }
  
  function(density, T) {
    DEBUG("v0: %f (micro Pa-sec), vEx: %f(micro Pa-sec)", vis0(T), visEx(density, T))
    (vis0(T) + visEx(density, T)) * 10^(-6)
  } 
})()


methaneProperties <- function(P, T) {
  Z <- methaneZ(P, T)
  density <- P * methaneM / (R * T) / methaneZ(P, T)
  viscosity <- methaneViscosity(density, T)
  list(density=density, viscosity=viscosity, Z=Z)
}

methaneZstd <- methaneZ(Pstp, Tstp)

