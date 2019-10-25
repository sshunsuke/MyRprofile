source("defaultSetting.R")

# =============================================================================
# Mukherjee & Brill method to predict pressure drop and flow regime in inclined
# two-Phase flow.
#
# The code in this file was written in accordance with the book titled 
# "Multiphase Flow in Wells" written by Brill & Mukherjee in 1999 although the 
# method was first published in Mukherjee & Brill (1985). 
# 
#   Mukherjee, H., Brill, J.P., 1985. Pressure drop correlations for inclined
#   two-phase flow. J. Energy Resour. Technol. Trans. ASME 107, 549–554.
# =============================================================================

# 
# densityMixS: mixture density with a consideration of slip

MukherjeeBrill <- list()

# DLNs
# flow regime
# holdup
# dPdL

MukherjeeBrill$calculate <- function(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle) {
	DLNs <- MukherjeeBrill$calculateDLNs(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle)
	flowRegime <- MukherjeeBrill$checkFlowRegime(DLNs)
	holdup <- MukherjeeBrill$holdup(DLNs, flowRegime)
	dPdL <- MukherjeeBrill$dPdL(DLNs, flowRegime, holdup)
	
}



# -----------------------------------------------------------------------------
# Calculate the dimensionless groups proposed by Duns & Ros. 
#
# Note:
#   densityG and viscosityG are not used in this function, but these are 
#   required to calculate dPdL. 
# -----------------------------------------------------------------------------
MukherjeeBrill$calculateDLNs = function(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, angle) {
	NLv <- vsL * (densityL / g / surfaceTension)^(0.25)
	NGv <- vsG * (densityL / g / surfaceTension)^(0.25)
	Nd <- D * sqrt(densityL * g / surfaceTension)
	NL <- viscosityL * (g / densityL / surfaceTension^3)^(0.25)
	
	# Gas velocity number for Slug/(Annular Mist) transition (4.130)
	NGvSM <- 10 ^ (1.401 - 2.694 * NL + 0.521 * NLv ^ 0.329)
	
	# Upflow: Bubble/Slug transition (4.128)
	x <- log10(NGv) + 0.940 + (0.074 * sin(angle)) - (0.855 * sin(angle)^2) + (3.695 * NL)
	NLvBS_up <- 10^x 
	
	# Downflow: Bubble/Slug transition (4.131)
	y <- 0.431 - (3.003 * NL) - (1.138 * log10(NLv) * sin(angle)) - (0.429 * log10(NLv)^2 * sin(angle)) + (1.132 * sin(angle))
	NGvBS <- 10^y
	
	# Downflow: Stratified (4.133)
	z <- 0.321 - ((0.017 * NGv) - 4.267 * sin(angle)) - (2.972 * NL) - (0.033 * log10(NGv)^2) - (3.925 * sin(angle)^2)
	NLvST <- 10^z
	
	data.frame(
		NLv = NLv,      # Liquid velocity number  (4.3)
		NGv = NGv,      # Gas velocity number     (4.4)
		Nd = Nd,        # Pipe diameter number    (4.5)
		NL = NL,        # Liquid viscosity number (4.6)
		
		NGvSM = NGvSM,            # Slug/(Annular Mist) transition (4.130)
		NGvBS = NGvBS,            # Downflow: Bubble/Slug transition (4.131)
		NLvBS_up = NLvBS_up,      # Upflow: Bubble/Slug transition (4.128)
		NLvST = NLvST,            # Downflow: Stratified (4.133)
		
		# Input values (these are used in the calculations of holdup and dPdL).
		vsG = vsG,
		vsL = vsL,
		D = D,
		densityG = densityG,
		densityL = densityL,
		viscosityG = viscosityG,
		viscosityL = viscosityL,
		surfaceTension = surfaceTension, 
		angle = angle   # [rad]
	)
}


# -----------------------------------------------------------------------------
# Determine flow regime.
# There are four flow regime types are defined in this function. 
#   1: Stratified
#   2: Annular
#   3: Slug
#   4: Bubbly
# -----------------------------------------------------------------------------
MukherjeeBrill$checkFlowRegime <- function(DLNs) {
  
  #fun <- function(DLNs) {
  fun <- function(NGv, NLv, angle, NGvSM, NGvBS, NLvBS_up, NLvST) {
    flowRegime <- 2  # annular 
    
    if (NGv > NGvSM) {
      return (flowRegime)  # Annular
    }
    
    if (angle > 0) {
      # Upfill
      if (NLv > NLvBS_up) {
        flowRegime <- 4  # bubbly
      } else {
        flowRegime <- 3  # slug
      }
    } else if (abs(angle) > deg2rad(-30)) {
      # Downhill
      if (NGv > NGvBS) {
        if (NLv > NLvST) {
          flowRegime <- 3  # Slug
        } else {
          flowRegime <- 1  # Stratified
        }
      } else {
        flowRegime <- 4    # bubbly
      }
    } else {
      # DownStratified
      if (NLv > NLvST) {
        if (NGv > NGvBS) {
          flowRegime <- 3    # Slug
        } else {
          flowRegime <- 4    # bubbly
        }
      } else {
        flowRegime <- 1  # Stratified
      }
    }
    
    return (flowRegime)
  }
  
  mapply(fun, DLNs$NGv, DLNs$NLv, DLNs$angle, DLNs$NGvSM, DLNs$NGvBS, DLNs$NLvBS_up, DLNs$NLvST)
}

# -----------------------------------------------------------------------------
# Estimate liquid holdup.
# -----------------------------------------------------------------------------

MukherjeeBrill$holdup <- function(DLNs, flowRegime) {
	co <- as.matrix( MukherjeeBrill$selectCoefficients(DLNs, flowRegime) )
	
	t1 <- co[1,] + (co[2,] * sin(DLNs$angle)) + (co[3,] * sin(DLNs$angle)^2) + (co[4,] * DLNs$NL^2)
	t2 <- DLNs$NGv^co[5,] / DLNs$NLv^co[6,]
	exp(t1 * t2)
}

MukherjeeBrill$selectCoefficients <- function(DLNs, flowRegime) {
	index <- ifelse(DLNs$angle > 0,
									1,
									ifelse((flowRegime == 1), 2, 3))
	MukherjeeBrill$coefficients[,index]
}

MukherjeeBrill$coefficients <- cbind(
	c(-0.380113, 0.129875, -0.119788,  2.343227, 0.475686, 0.288657),   # Up
	c(-1.330282, 4.808139,  4.171584, 56.262268, 0.079951, 0.504887),   # DownStratified
	c(-0.516644, 0.789805,  0.551627, 15.519214, 0.371771, 0.393952)    # DownOther
)
colnames(MukherjeeBrill$coefficients) <- c("Up", "DownStratified", "DownOther")
rownames(MukherjeeBrill$coefficients) <- paste("C", 1:6, sep="")
	

# -----------------------------------------------------------------------------
# dPdL
# -----------------------------------------------------------------------------

# Usage: MukherjeeBrill$ffRatio(HR)
MukherjeeBrill$ffRatio <- (function(){
	# Table 4.5
	correspondence <- cbind(
		c(1.00, 0.98, 1.20, 1.25, 1.30, 1.25, 1.00,  1.00),    # fR
		c(0.01, 0.20, 0.30, 0.40, 0.50, 0.70, 1.00, 10.00)     # HR
	)
	colnames(correspondence) <- c("fR", "HR")
	approxfun(correspondence[,"HR"], correspondence[,"fR"])
})()


MukherjeeBrill$dPdL <- function(DLNs, flowRegime, HL) {
	dPdL <- NA

	if (flowRegime == 1) {
		# Stratified
		
		# delta: angle related to liquid level (shown in Fig 4.20)
		fun <- function(delta) { 1/(2*pi) * (delta - sin(delta)) - HL }  # (4.147)
		f <- uniroot(fun, c(0,2*pi))
		delta <- f$root
		
		# Flow area of each phase
		AL <- circle$area(DLNs$D / 2) * HL         # (4.147)
		AG <- circle$area(DLNs$D / 2) * (1 - HL)
		
		
		# 4.146 for hL/d
		#   4.150 and 4.151
		dhG <- DLNs$D * (2*pi - (delta - sin(delta))) / (2*pi - delta + 2 * sin(delta/2))    # (4.150)
		dhL <- DLNs$D * (delta - sin(delta)) / (delta + 2 * sin(delta / 2))                  # (4.151)
		
		# Perimeter
		P <- DLNs$D * pi
		PG <- (1 - delta / (2*pi)) * P    # 4.149
		PL <- P - PG                      # 4.148
		
		# Shear stress
		vG <- DLNs$vsG / (1-HL)           # 4.157
		vL <- DLNs$vsL / HL               # 4.156
		
		ReG <- DN$Reynolds(DLNs$densityG, vG, dhG, DLNs$viscosityG)   # (4.155)
		ReL <- DN$Reynolds(DLNs$densityL, vL, dhL, DLNs$viscosityL)   # (4.154)
		
		fDG <- FCP$fD.Blasius(ReG)
		fDL <- FCP$fD.Blasius(ReL)
		
		shearStressG <- fDG * DLNs$densityG * vG^2 / (2*g)    # 4.153
		shearStressL <- fDL * DLNs$densityL * vL^2 / (2*g)    # 4.152
		
		# dPdL (4.144)
		dPdL <- - (shearStressL * PL + shearStressG * PG) - (DLNs$viscosityL * AL + DLNs$viscosityG * AG) * g * sin(DLNs$angle)
	} else {
		vmix <- DLNs$vsG + DLNs$vsL
		HLnoslip <- DLNs$vsL / vmix
		
		densityMixS <- DLNs$densityG * (1 - HL) + DLNs$densityL * HL                    # (3.22)
		densityMixN <- DLNs$densityG * (1 - HLnoslip) + DLNs$densityL * HLnoslip        # (3.23)
		viscosityMixS <- DLNs$viscosityG * (1 - HL) + DLNs$viscosityL * HL              # (3.19)
		viscosityMixN <- DLNs$viscosityG * (1 - HLnoslip) + DLNs$viscosityL * HLnoslip  # (3.21)
		
		ReN <- DN$Reynolds(densityMixN, vmix, DLNs$D, viscosityMixN)

		if (flowRegime == 2) {
			# Annular
			fn <- FCP$fD.Blasius(ReN)             # no-slip friction factor
			HR <- HLnoslip / HL                   # (4.140)
			fR <- MukherjeeBrill$ffRatio(HR)      # friction factor ratio
			fD <- fn * fR                         # (4.141)
			Ek <- densityMixS * vmix * vsG / p    # (4.53)  (4.137)
			dPdL <- (fD * densityMixS * vmix^2 / (2 * D) + densityMixS * g * sin(angle)) / (1 - Ek)  # (4.139)
		} else if (flowRegime == 3 | flowRegime == 4) {
			# Slug or Bubble
			fD <- FCP$fD.Blasius(ReN)
			Ek <- densityMixS * vmix * DLNs$vsG / DLNs$  # (4.137)
			dPdL <- (fD * densityMixS * vmix^2 / (2 * D) + densityMixS * g * sin(angle)) / (1 - Ek)  # (4.136)
		} else {
			# error!
		}
	}
	
	dPdL
}


MukherjeeBrill$dPdLSlugBubble <- function(){
	Re <- densityMixN * vmix * DLNs$D / viscosityMixN
	fD <- FCP$fD.Blasius(Re)
	Ek <- densityMixS * vmix * vsG / p  # (4.137)
	dPdL <- (fD * densityMixS * vmix^2 / (2 * D) + densityMixS * g * sin(angle)) / (1 - Ek)  # (4.136)
}





# =============================================================================

# =============================================================================


isAnnularMist = function(DLNs) { DLNs$NGv > DLNs$NGvSM }

MukherjeeBrill$coefficients


# vsL, vsG, D, densityL, surfaceTension, viscosityL, angle
#
# angle (rad)


# Example 3.2
#   d = 6", vsG = 3.86 ft/sec, vsL = 3.97 ft/sec
# Example 4.1
#   densityG = 5.88 lbm/ft3, densityL = 47.61 lbm/ft3, 
# Example 4.8
#   viscosityG = 0.016 cp, viscosityL = 0.97 cp, surfaceTension = 8.41 dynes/cm, angle = 90 deg
ratio <- lbm2kg(1) / (ft2m(1)^3)
exam <- MukherjeeBrill$calculateDLNs(ft2m(3.86), ft2m(3.97), inch2m(6), 5.88*ratio, 47.61*ratio, UC$cP2Pas(0.016), UC$cP2Pas(0.97), UC$dynpcm2Npm(8.41), pi/2)
fr <- MukherjeeBrill$checkFlowRegime(exam)
hol <- MukherjeeBrill$holdup(exam, fr)
MukherjeeBrill$dPdL(exam, fr, hol)

exam2 <- MukherjeeBrill$calculateDLNs(ft2m(c(3.86, 0.1)), ft2m(3.97), inch2m(6), 5.88*ratio, 47.61*ratio, UC$cP2Pas(0.016), UC$cP2Pas(0.97), UC$dynpcm2Npm(8.41), c(pi/2, -0.1))
fr2 <- MukherjeeBrill$checkFlowRegime(exam2)
MukherjeeBrill$selectCoefficients(exam2, fr2)
MukherjeeBrill$holdup(exam2, fr2)













if (FALSE) {
	exam <- MukherjeeBrill$calculateDLNs(3.86, 3.97, 0.5, 5.88, 47.61, 0.016, 0.97, 8.41, pi/2)
	
	a <- MukherjeeBrill$calculateDLNs(10,1,inch2m(3), 1000, 0.07, 0.001026800601060196, 0)
	MukherjeeBrill$checkFlowRegime(a)
	
	# Dimensionless Groups proposed by Duns & Ros
	
	# Liquid velocity number
	NLv = vsL * (densityL / g / surfaceTension)^(0.25)
	
	# Gas velocity number
	NGv = vsG * (densityL / g / surfaceTension)^(0.25)
	
	# Pipe diameter number
	Nd = D * sqrt(densityL * g / surfaceTension)
	
	# Liquid viscosity number
	NL = viscosityL * (g / densityL / surfaceTension^3)^(0.25)
	
	
	
	
	
	x = log10(NGv) + 0.940 + 0.074 * sin(angle) - 0.855 * sin(angle)^2 + 3.695 * NL
	
	NL
	
	
	{
		if (DLNs$angle <= 0) {
			if (flowRegime == 1) {
				return(MukherjeeBrill$coefficients[,"DownStratified"])
			} else {
				return(MukherjeeBrill$coefficients[,"DownOther"])
			}
		}
		MukherjeeBrill$coefficients[,"Up"]
	}
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

