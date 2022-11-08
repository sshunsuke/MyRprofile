source("defaultSetting.R")


xyTableFromFunction <- function(fun, x, ...) {
	y <- fun(x, ...)
	cbind(x,y)
}



hydraulicDiameter <- function(A, S) {
  4 * A / S
}




# Drift Flux Model. ----
#
#   1. vmix = vsg + vsl
#   2. vg = C0 * vmix + vd
#   3. voidRatio = vsg / vg
#

cpf.driftFluxModel <- function(C0, Vd, vsg, vsl) {
  vmix <- vsg + vsl
  vg <- C0 * vmix + Vd
  voidRatio <- vsg / vg
  
  Hl <- 1 - voidRatio
  vl <- vsg / Hl
  
  c(vg=vg, vl=vl, Hl=Hl)
}


cpf.dfm.C0 <- list(
  bubble = function(densityG, densityL) {
    
  }
  
  
)


# ----------
# Bahadori

A1 = c( -1.619063838,           -1.4673416207,         -1.3064287345,         -1.5068329304         )
B1 = c( -6.0440641629*10^(-2),  -9.4579004057*10^(-3),  3.7689223988*10^(-2),  1.6368010607*10^(-1) )
C1 = c(  1.2992412636*10^(-3),  -9.0991682769*10^(-4), -5.7653162538*10^(-3), -1.9272634157*10^(-2) )
D1 = c( -1.0480516067*10^(-5),   2.1036093111*10^(-5),  1.5920740612*10^(-4),  6.081932874*10^(-4)  )

A2 = c( -5.675424778 *10^(-2),   5.6420129717*10^(-3),  2.3856855469*10^(-2),  2.1494907991*10^(-1) )
B2 = c(  1.1266206576*10^(-3),  -8.1324216389*10^(-3), -1.7666292295*10^(-2), -7.236908659 *10^(-2) )
C2 = c( -4.5476251244*10^(-5),  3.6013086233 *10^(-4),  1.2527015231*10^(-3),  6.2521396767*10^(-3) )
D2 = c(  5.4011484658*10^(-7),  -5.0991959691*10^(-6), -2.9095395202*10^(-5), -1.7702985302*10^(-4) )

A3 = c(  1.287145175 *10^(-3),  -1.0914548287*10^(-3), -2.2442814135*10^(-3), -1.2571952312*10^(-2) )
B3 = c( -3.1321987972*10^(-5),   3.1336503528*10^(-4),  7.879778678 *10^(-4),  3.6872340957*10^(-3) )
C3 = c(  1.3585299744*10^(-6),  -1.3704000608*10^(-5), -5.6411819564*10^(-5), -3.1830356934*10^(-4) )
D3 = c( -1.6951529528*10^(-8),   1.9239332368*10^(-7),  1.3264010427*10^(-6),  9.0141402462*10^(-6) )

A4 = c( -1.1238180847*10^(-5),   1.4969220576*10^(-5),  3.3668468005*10^(-5),  1.4981731667*10^(-4) )
B4 = c(  3.574027836 *10^(-7),  -3.3674152566*10^(-6), -9.9103300792*10^(-6), -4.2562432156*10^(-5) )
C4 = c( -1.5140460085*10^(-8),   1.4658761606*10^(-7),  7.1509244179*10^(-7),  3.6706099161*10^(-6) )
D4 = c(  1.8649575376*10^(-10), -2.0517177043*10^(-9), -1.6973934089*10^(-8), -1.0392224103*10^(-7) ) 


ce_100C = cbind(c(A1[1], B1[1], C1[1], D1[1]),
                c(A2[1], B2[1], C2[1], D2[1]),
                c(A3[1], B3[1], C3[1], D3[1]),
                c(A4[1], B4[1], C4[1], D4[1]))

rownames(ce_100C) <- c('A', 'B', 'C', 'D')
colnames(ce_100C) <- c('alpha', 'beta', 'gamma', 'theta')


calc_thickness <- function(d, tc, ce) {
  alpha = ce['A','alpha'] + ce['B','alpha']/tc + ce['C','alpha'] / tc^2 + ce['D','alpha'] / tc^3
  beta = ce['A','beta'] + ce['B','beta']/tc + ce['C','beta'] / tc^2 + ce['D','beta'] / tc^3
  gamma = ce['A','gamma'] + ce['B','gamma']/tc + ce['C','gamma'] / tc^2 + ce['D','gamma'] / tc^3
  theta = ce['A','theta'] + ce['B','theta']/tc + ce['C','theta'] / tc^2 + ce['D','theta'] / tc^3
  
  thickness = exp( alpha + beta / d + gamma / d^2 + theta / d^3)
  
  cat( sprintf('alpha = %e, beta = %e, gamma = %e, theta = %e \n', alpha, beta, gamma, theta) )
  cat( sprintf('thickness = %e \n', thickness) )
  thickness
}

calc_thickness(0.25, 0.04, ce_100C)


# ----------

  
# dotQ = dotm * cp * (Tout - Tin) = LMTD * htc * pi * D * dL 

# dotm * cp = htc * (pi * D * dL) * 1 / (log((Tamb - Tin) / (Tamb - Tout)))


# log((Tamb - Tin) / (Tamb - Tout)) = htc * (pi * D * dL) / dotm / cp
# (Tamb - Tin) / (Tamb - Tout) = exp( htc * (pi * D * dL) / dotm / cp )
# Tamb - Tout =  (Tamb - Tin)  /  (exp( htc * (pi * D * dL) / dotm / cp ))


f <- HeatExchange$Tprofile$funConstantAmbientT(20, 4, 450, UC$inch2m(14), 94, 2000)
plot(f, to=1000)


xyTableFromFunction(f, 1:100*10)



# -----------------------------------------------------------------------------