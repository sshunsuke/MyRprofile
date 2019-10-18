# =============================================================================
# Mukherjee & Brill.
# =============================================================================


MukherjeeBrill <- list()

MukherjeeBrill$coefficients <- cbind(
	c(-0.380113, 0.129875, -0.119788,  2.343227, 0.475686, 0.288657),
	c(-1.330282, 4.808139,  4.171584, 56.262268, 0.079951, 0.504887),
	c(-0.516644, 0.789805,  0.551627, 15.519214, 0.371771, 0.393952)
)
colnames(MukherjeeBrill$coefficients) <- c("Up", "DownStratified", "Down")
rownames(MukherjeeBrill$coefficients) <- paste("C", 1:6, sep="")

# calculate the dimensionless groups proposed by Duns & Ros
MukherjeeBrill$calculateDLNs = function(vsG, vsL, D, densityL, surfaceTension, viscosityL, angle) {
	NLv <- vsL * (densityL / g / surfaceTension)^(0.25)
	NGv <- vsG * (densityL / g / surfaceTension)^(0.25)
	Nd <- D * sqrt(densityL * g / surfaceTension)
	NL <- viscosityL * (g / densityL / surfaceTension^3)^(0.25)
	
	# Gas velocity number for Slug/(Annular Mist) transition (4.130)
	NGvSM <- 10 ^ (1.401 - 2.694 * NL + 0.521 * NLv ^ 0.329)
	
	# Upflow: Bubble/Slug transition (4.128)
	x <- log(NGv) + 0.940 + 0.074 * sin(angle) - 0.855 * sin(angle)^2 + 3.695 * NL
	NLvBS_up <- 10^x 
	
	# Downflow: Bubble/Slug transition (4.131)
	y <- 0.431 - (3.003 * NL) - (1.138 * log(NLv) * sin(angle)) - (0.429 * log(NLv)^2 * sin(angle)) + (1.132 * sin(angle))
	NGvBS <- 10^y
	
	# Downflow: Stratified (4.133)
	z <- 0.321 - ((0.017 * NGv) - 4.267 * sin(angle)) - (2.972 * NL) - (0.033 * log(NGv)^2) - (3.925 * sin(angle)^2)
	NLvST <- 10^z
	
	list(
		NLv = NLv,      # Liquid velocity number  (4.3)
		NGv = NGv,      # Gas velocity number     (4.4)
		Nd = Nd,        # Pipe diameter number    (4.5)
		NL = NL,        # Liquid viscosity number (4.6)
		
		NGvSM = NGvSM,            # Slug/(Annular Mist) transition (4.130)
		NLvBS_up = NLvBS_up,      # Upflow: Bubble/Slug transition (4.128)
		NGvBS = NGvBS,  # Downflow: Bubble/Slug transition (4.131)
		NLvST = NLvST,  # Downflow: Stratified (4.133)
		
		angle = angle   # [rad]
	)
}

# Check flow regime.
# There are four flow regime types are defined in this function. 
#   1: Stratified
#   2: Annular
#   3: Slug
#   4: Bubbly
MukherjeeBrill$checkFlowRegime <- function(DLNs) {
	flowRegime <- 2  # annular 
	
	if (DLNs$NGv > DLNs$NGvSM) {
		return (flowRegime)
	}
	
	if (DLNs$angle > 0) {
		# Upfill
		if (DLNs$NLv > DLNs$NLvBS_up) {
			flowRegime <- flowRegime <- 4  # bubbly
		} else {
			flowRegime <- flowRegime <- 3  # slug
		}
	} else if (abs(DLNs$angle) > deg2rad(-30)) {
		# Downhill
		if (DLNs$NGv > DLNs$NGvBS) {
			if (DLNs$NLv > DLNs$NLvST) {
				flowRegime <- flowRegime <- 3  # Slug
			} else {
				flowRegime <- flowRegime <- 1  # Stratified
			}
		} else {
			flowRegime <- flowRegime <- 4    # bubbly
		}
	} else {
		# DownStratified
		if (DLNs$NLv > DLNs$NLvST) {
			if (DLNs$NGv > NGvBS) {
				flowRegime <- flowRegime <- 3  # Slug
			} else {
				flowRegime <- flowRegime <- 4    # bubbly
			}
		} else {
			flowRegime <- flowRegime <- 1  # Stratified
		}
		
	}
	
	return (flowRegime)
}


isAnnularMist = function(DLNs) { DLNs$NGv > DLNs$NGvSM }

MukherjeeBrill$coefficients


# vsL, vsG, D, densityL, surfaceTension, viscosityL, angle
#
# angle (rad)

a <- MukherjeeBrill$calculateDLNs(10,1,inch2m(3), 1000, 0.07, 0.001026800601060196, 0)
MukherjeeBrill$checkFlowRegime(a)








if (FALSE) {
	# Dimensionless Groups proposed by Duns & Ros
	
	# Liquid velocity number
	NLv = vsL * (densityL / g / surfaceTension)^(0.25)
	
	# Gas velocity number
	NGv = vsG * (densityL / g / surfaceTension)^(0.25)
	
	# Pipe diameter number
	Nd = D * sqrt(densityL * g / surfaceTension)
	
	# Liquid viscosity number
	NL = viscosityL * (g / densityL / surfaceTension^3)^(0.25)
	
	
	
	
	
	x = log(NGv) + 0.940 + 0.074 * sin(angle) - 0.855 * sin(angle)^2 + 3.695 * NL
	
	NL
	
}





# =============================================================================

xyTableFromFunction <- function(fun, x, ...) {
	y <- fun(x, ...)
	cbind(x,y)
}



hydraulicDiameter <- function(A, S) {
  4 * A / S
}




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


f <- HeatExchange$Tprofile$funConstantAmbientT(20, 4, 450, UC$inch2m(14), 94, 2000)
plot(f, to=1000)


xyTableFromFunction(f, 1:100*10)

