

xyTableFromFunction <- function(fun, x, ...) {
	y <- fun(x, ...)
	cbind(x,y)
}



hydraulicDiameter <- function(A, S) {
  4 * A / S
}



# =============================================================================
# Friction Factor.
# =============================================================================

# Fanning Friction Factor
cpf.ff <- (function(){
  
  colebrook_ <- function(roughness, D, Re, interval=c(0, 5)) {
    core_ <- function(roughness, D, Re, interval, warn=TRUE) {
      if (Re <= 4000 && (cpf.ff$WARNING == TRUE)) { warning("Re <= 4000 !") }
      
      fun <- function(f) {
        (1 / sqrt(f)) + 4 * log10( roughness / D / 3.71 + 1.256 / Re / sqrt(f))
      }
      f <- uniroot(fun, interval)  # newton.method
      f$root
    }
    
    wrap_ <- function(roughness_, D_, Re_){
      core_(roughness_, D_, Re_, interval)
    }
    mapply(wrap_, roughness, D, Re)
  }
  
  list(
    WARNING = TRUE,
    
    Colebrook = colebrook_
  )
})()



# Borda-Carnot's formula
#   dE = lossCoefficient * density * (vout - vin)^2 / 2
#     v: flow velocity
#     density: 
#     Din: diameter before expansion
#     Dout: diameter after expansion
#     eta: pressure recovery factor (optional parameter)
BordaCarnot = function(vin, density, Din, Dout, eta){
	if (missing(eta)) { eta <- 0 }
	lossCoefficient = (1 - eta) * (1 - (Din / Dout)^2)^2
	lossCoefficient * density * vin^2 / 2
}

# BordaCarnotHead



# =============================================================================
# Drift Flux Model.
#
#   1. vmix = vsg + vsl
#   2. vg = C0 * vmix + vd
#   3. voidRatio = vsg / vg
#
# =============================================================================

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






# LMTD (logarithmic mean temperature difference)
HeatExchange <- list(
  
  
  # Heat Transfer Rate [W = J/sec]
  heatTransferRate = function(LMTD, heatTransferCoefficient, SurfaceArea) {
    LMTD * heatTransferCoefficient * SurfaceArea
  },
  
  Tprofile = list(
    funConstantAmbientT = function(Tin, Tamb, htc, D, mdot, cp) {
      f <- function(x) {
        exp(- htc * (pi * D) / (mdot * cp) * x) * (Tin - Tamb) + Tamb
      }
      return(f)
    }
  ),
  

  LMTD = list(
    
    # Meaning of the symbols:
    #   dT   in         out
    #
    #   T1   in ------> out
    #   T2  out <------ in
    counterCurrent = function(T1_in, T1_out, T2_in, T2_out) {
      # n <- (T1_in - T2_out) - (T1_out - T2_in)
      # d <- log((Tain - Tbout) / (Taout - Tbin))
      # n / d
      
      dT_in  <- T1_in - T2_out
      dT_out <- T1_out - T2_in
      
      (dT_in - dT_out) / log(dT_in / dT_out)
    },
    
    
    # Positive value means heat flow from line-1 to line-2. 
    # Nagative value means heat flow from line-2 to line-1. 
    
    # Meaning of the symbols:
    #   dT   in         out
    #
    #   T1   in ------> out
    #   T2   in ------> out
    coCurrent = function(T1_in, T1_out, T2_in, T2_out) {
      if (T1_in > T2_in) {
        if (T1_out < T2_out) {
          stop("'T1_out < T2_out' despite 'T1_in > T2_in'. ")
        } else if (T1_out > T1_in) {
          stop("'T1_out > T1_in' despite 'T1_in > T2_in'. ")
        } else if (T2_out < T2_in) {
          stop("'T2_out < T2_in' despite 'T1_in > T2_in'. ")
        }
      } else if (T1_in < T2_in) {
        if (T1_out > T2_out) {
          stop("'T1_out > T2_out' despite 'T1_in < T2_in'. ")
        } else if (T1_out < T1_in) {
          stop("'T1_out < T1_in' despite 'T1_in < T2_in'. ")
        } else if (T2_out > T2_in) {
          stop("'T2_out > T2_in' despite 'T1_in < T2_in'. ")
        }
      } else {
        if (T1_out != T2_out) {
          stop("'T1_out != T2_out' despite 'T1_in == T2_in'. ")
        }
      }

      dT_in  <- T1_in - T2_in
      dT_out  <- T1_out - T2_out
      (dT_in - dT_out) / log(dT_in / dT_out)
    },
    
    # LMTD for a pipe flow with constant ambient T
    #         Tamb
    # T   in ------> out
    constantAmbientT = function(Tin, Tout, Tamb) {
      HeatExchange$LMTD$coCurrent(Tin, Tout, Tamb, Tamb)
    }
  )
  
)



  
# dotQ = dotm * cp * (Tout - Tin) = LMTD * htc * pi * D * dL 

# dotm * cp = htc * (pi * D * dL) * 1 / (log((Tamb - Tin) / (Tamb - Tout)))


# log((Tamb - Tin) / (Tamb - Tout)) = htc * (pi * D * dL) / dotm / cp
# (Tamb - Tin) / (Tamb - Tout) = exp( htc * (pi * D * dL) / dotm / cp )
# Tamb - Tout =  (Tamb - Tin)  /  (exp( htc * (pi * D * dL) / dotm / cp ))


f <- HeatExchange$Tprofile$funConstantAmbientT(20, 4, 450, UC$inch2m(14), 18, 2000)
plot(f, to=1000)


xyTableFromFunction(f, 1:100*10)

