# Comment Stype: https://ourcodingclub.github.io/tutorials/etiquette/

# Nomenclatures for properties:
#
# fD: Darcy friction factor
# g: gravity Acceleration [m/s2]
# r: radius 
# v: velocity [m/sec]
# x: molar fraction
#
# A: area [m2]
# D: diameter [m]
# L: length [m]
# P: pressure [Pa]
# R: Gas Constant [J/K-mol]
# T: temperature [K]
# V: volume [m3]
# 
# 

# _____________________________________________________________________________


# *****************************************************************************
# ** UTIL (Utilities) ** ----
# *****************************************************************************

if (exists('SMD') == TRUE) { detach(UTIL) }

UTIL <- list(
  
  # ***************************
  # * Constant values ----
  # ***************************
  
  g = 9.8,                      # Gravity Acceleration (m/s2)
  R = 8.3144621,                # Gas Constant (J/K-mol)
  # R <- 8.314471               # Moldover et al. (1988)
  
  # Na <- 6.022140857 * 10^23     # Avogadro constant
  kB = 1.38064852 * 10^(-23),   # Boltzmann constant (J/K)
  
  # STP (Standard Temperature and Pressure) - old 
  Pstp = 100000,    # (Pa)
  Tstp = 273.15,    # (K)
  
  # SATP (Standard Ambient Temperature and Pressure)
  Psatp = 101325,    # (Pa)
  Tsatp = 298.15,    # (K)
  
  # NTP (Normal Temperature and Pressure)
  Pntp = 101325,     # (Pa)
  Tntp = 293.15,     # (K)
  
  
  # ***********************************
  # * Properties of simple shapes ----
  # ***********************************
  
  circle = list(
    area = function(r) { r * r * pi },
    circumference = function(r) { 2 * r * pi }
  ),
  
  sphere = list(
    volume      = function(r) { (4 * pi * r^3) / 3 },
    surfaceArea = function(r) { 4 * pi * r^2 }
  ),
  
  elipsoid = list(
    volume      = function(a,b,c) { 4/3 * pi * a * b * c }
  ),
  
  # Sauter Mean Diameter (vD: Vector of diameters)
  SMD = function(vD) { sum(vD^3) / sum(vD^2) },
  
  shape2d = list(
    extent = function(A, Arect) { A / Arect },
    solidity = function(A, Ahull) { area / Ahull },
    circularity = function(A, peri) { 4 * pi * A / (peri^2) }
  ),
  
  
  # ***********************************
  # * IO Clipboard ----
  # ***********************************
  
  rcbmat = function(header, ...) {
    if (missing(header)) { header = FALSE }
    utils::read.table(file="clipboard", header=header, ...)
  },
  
  wcbmat = function(data, header, sep="\t", size=128, row.names=FALSE, qmethod="double",
                    col.names=ifelse(header && row.names, NA, header), ...) {
    if (missing(header)) { header = FALSE }
    
    if (typeof(size) != "double" && typeof(size) != "integer") {
      stop("type of 'size' must be 'double' or 'integer'.")
    }
    fn = sprintf("clipboard-%d", size)
    
    utils::write.table(data, file=fn, sep=sep, row.names=row.names,
                       col.names=col.names, qmethod=qmethod, ...)
  },
  
  
  # ***********************************
  # * Data ----
  # ***********************************
  
  # Vector
  v.indexClosestValue = function(vec, x) {
    diff <- abs(vec-x)
    which( diff == min(diff) )
  },
  
  v.closestValue = function(vec, x) {
    vec[ v.indexClosestValue(vec, x) ]
  },
  
  # ***************
  # Matrix
  # ***************
  m.crev = function(mat) { mat[,ncol(mat):1] },
  m.rrev = function(mat) { mat[nrow(mat):1,] },
  
  # ***************
  # * Data Frame
  # ***************
  df.orderBy = function(df, colname, decreasing=FALSE) {
    if (missing(colname)) { stop("'colname' is not specified.") }
    df[order(df[,colname], decreasing=decreasing),]
  },

  # ***********************************
  # * Probability
  # ***********************************
  
  # Create a matrix of information of cumulative distribution.
  # 
  # If you want to create a graph of cumulative distribution, you should use
  # ecdf() function.
  #   > fCP <- ecdf( c(3,76,58,24,100,1) ); plot(fCP)
  cum.probability = function(values, decreasing=FALSE) {
  	len <- length(values)
  	cbind(X = sort(values, decreasing=decreasing), cum.prob = (1:len)/len)
  },
  
  cum.probability.guessX = function(cp_mat, cp) {
    diff_cp = abs(cp_mat[,"cum.prob"] - cp)
    
    rn = which(diff_cp == min(diff_cp))
    if (min(diff_cp) == 0) {
      return( as.numeric(cp_mat[rn,"X"]) )
    }
    
    if (length(rn) == 1) {
      rn2 = ifelse(cp_mat[rn,"cum.prob"] > cp, rn-1, rn+1)
    } else if (length(rn) == 2) {
      rn2 = rn[2]
      rn = rn[1]
    } else {
      stop()
    }
    
    p = sort(cp_mat[rn:rn2,"cum.prob"])
    x = sort(cp_mat[rn:rn2,"X"])
    
    x[1] + (x[2] - x[1]) * (cp - p[1]) / (p[2] - p[1])
  },
  
  # ***********************************
  # * Root-finding algorithm ----
  # ***********************************
  
  # Newton-Raphson method
  newtonRaphson = function(fun, dFun, fD_0, tol=1e-10, itMax=10) {
    it <- 0
    fD_n <- fD_0
    d <- fun(fD_n) / dFun(fD_n)
    
    while (abs(d) >= tol) {
      it <- it + 1
      if (it > itMax) {
        stop("Calculation did not converge.")
      }
      
      d <- fun(fD_n) / dFun(fD_n)
      fD_n <- fD_n - d
    }
    
    fD_n
  },
  
  
  # Bisection Method
  bisection = function(f, rangeFrom, rangeTo, itMax = 100, tol = 1e-7) {
    # Check arguments.
    if (rangeFrom >= rangeTo) {
      stop("'rangeFrom' must be smaller than 'rangeTo'.")
    } else if (f(rangeFrom) * f(rangeTo) >= 0) {
      stop('the sign of f(rangeFrom) must be different from that of f(rangeTo)')
    } 
    
    a <- rangeFrom
    b <- rangeTo
    c <- a
    
    for (i in 1:itMax) {
      c <- (a + b) / 2 # Calculate midpoint
      
      # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the 
      # function and return the root.
      if (abs(f(c)) < tol) {
        break
      }
      
      # If another iteration is required, 
      # check the signs of the function at the points c and a and reassign
      # a or b accordingly as the midpoint to be used in the next iteration.
      if (sign(f(c)) == sign(f(a))) {
        a <- c
      } else {
        b <- c
      }
    }
    
    if (i >= itMax) {
      print('Too many iterations')
    }
    
    return(c)
  }
  
  
  
)

if (exists('SMD') == FALSE) { attach(UTIL) }



# *******************************
# ** UC (Unit Conversion) ** ----
# *******************************
if (exists('inch2m') == TRUE) { detach(UC) }

UC <- list(
  # Angle
  rad2deg = function(rad) { rad * 180 / pi },
  deg2rad = function(deg) { deg * pi / 180 },
  
  # Length
  inch2m = function(inch) { inch * 0.0254 },
  m2inch = function(m) { m / 0.0254 },
  ft2m   = function(ft) { ft * 0.3048 },
  m2ft   = function(m) { m / 0.3048 },
  
  # Volume
  bbl2m3 = function(bbl)  { bbl * 0.158987294928 },   # (oil) barrel -> m3
  m32bbl = function(m3)   { m3 / 0.158987294928},     # m3 -> (oil) barrel
  bbl2gal = function(bbl) { bbl * 42 },               # (oil) barrel -> (US) bbl
  gal2bbl = function(gal) { gal / 42 },               # (US) bbl -> (oil) barrel
  
  # Temperature
  K2C = function(K) { K - 273.15 },
  C2K = function(C) { C + 273.15 },
  F2C = function(F) { (F - 32) * 5 / 9 },
  C2F = function(C) { (C * 9/5) + 32 },
  F2R = function(F) { F + 459.67 },
  R2F = function(R) { R - 459.67 },
  
  # Pressure
  psi2Pa = function(psi) { psi * 6894.76 },  # psi (pound-force per square inch) -> Pa
  Pa2psi = function(Pa)  { Pa / 6894.76 },   # Pa -> psi
  atm2Pa = function(atm) { atm * 101325 },   # atm -> Pa
  Pa2atm = function(Pa)  { Pa / 101325 },    # Pa -> atm
  
  # Weight
  lbm2kg = function(lbm) { lbm * 0.45359237 },       # pound-mass -> kg
  kg2lbm = function(kg) { kg / 0.45359237 },         # kg -> pound-mass
  
  # Density
  kgcbm2lbcbft = function(kgcbm) { kgcbm * 0.062428},  # kg/m3 -> lb/ft3 (pound per cubic foot)
  
  
  # Angular velocity <-> Frequency
  radps2Hz = function(radps) { radps / (2 * pi) },  # rad/s -> Hz
  Hz2radps = function(Hz) { Hz * (2 * pi) },        # Hz -> rad/s
  
  # 
  
  # Viscosity
  cP2Pas = function(cP) {cP / 1000},
  Pas2cP = function(Pas) {Pas * 1000},
  
  # Surface Tension
  Npm2dynpcm = function(Npm) {Npm * 1000},
  dynpcm2Npm = function(dynpcm) {dynpcm / 1000}
)

if (exists('inch2m') == FALSE) { attach(UC) }





# ************************************
# ** DN (Dimensionless Number) ** ----
# ************************************
if (exists('Reynolds') == TRUE) { detach(DN) }

DN <- (function(){
  
  # Eotvos number (also called the Bond number)
  #   dDensity: difference in density f the two phase [kg/m3]
  Eo_Bo <- function(dDensity, L, surfaceTension) {
    (dDensity * g * L^2) / surfaceTension
  }
  
  
  list(
    
    # Bond number (Bo)
    Bond = Eo_Bo,
    
    # Eotvos number (Eo)
    Eotvos = Eo_Bo,
    
    # Mo
    Morton = function(cDensity, cViscosity, surfaceTension, dDensity) {
      if (missing(dDensity)) {
        dDensity <- cDensity
      }
      (g * cViscosity^4 * dDensity) / (cDensity^2 * surfaceTension^3) 
    },
    
    
    # Nu := "convective heat transfer" / "conductive heat transfer"
    #   heatTransferCoefficient: W/m2-K, thermalConductivity: W/m-K
    Nusselt = function(heatTransferCoefficient, thermalConductivity, L) {
    	heatTransferCoefficient * L / thermalConductivity
    },

    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Ratios between diffusions of momentum, heat or mass
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Pr := momentum / heat
    #   specificHeat: J/Kg-K, viscosity: N-s/m2, thermal conductivity: W/m-K
    Prandtl = function(specificHeat, viscosity, thermalConductivity) {
      (specificHeat * viscosity) / thermalConductivity
    },
    
    # Sc := momentum / mass
    Schmidt = function(viscosity, density, diffusionCoefficient) {
      viscosity / (density * diffusionCoefficient)
    },
    
    # Le := heat / mass
    Lewis = function(thermalConductivity, density, specificHeat, diffusionCoefficient) {
    	thermalConductivity / (density * specificHeat * diffusionCoefficient)
    },
    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Ratio of different forces
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    # Fr
    Froude = function(v, L) {
      v / (L * g)
    },
    
    # Gr := "bouyancy" / "viscous force"
    # 
    # Heat Transfer:
    #   beta = thermal expansion coefficient
    #   dT_dC = temperature difference
    # Mass Transfer:
    #   
    Grashof = function(L, beta, dT_dC, kinematicViscosity) {
      (g * beta * dT_dC * L^3) / (kinematicViscosity^2)
    },
    
    # We: [= "inertia" / "surface tension"]
    Weber = function(density, v, L, surfaceTension) {
      density * (v^2) * L / surfaceTension
    },
    
    # Re: [= "inertia" / "viscous force"]
    Reynolds = function(density, v, L, viscosity) {
      density * v * L / viscosity
    }
  )
})()

if (exists('Reynolds') == FALSE) { attach(DN) }







# **********************************************
# ** DEBUG_LOG (Functions for debug log) ** ----
# **********************************************
DEBUG_LOG <- (function() {
  FILE        <- "log/rlog.txt"
  
  LEVEL_FATAL <- 1
  LEVEL_ERROR <- 2
  LEVEL_WARN  <- 3
  LEVEL_INFO  <- 4
  LEVEL_DEBUG <- 5
  LEVEL_TRACE <- 6
  LEVEL       <- LEVEL_DEBUG
  CONSOLE <- TRUE
  
  INDENT <- (function(n) {
    s <- "  "
    l <- numeric(n)
    l[1] <- ""
    for (i in 2:n) {
      l[i] <- paste(l[i-1], s, sep="")
    }
    l
  })(10)
  
  OUTPUT <- function(level, indent, msg, console=TRUE) {
    s <- paste(indent, msg, sep="")
    if (console && DEBUG_LOG$CONSOLE) { cat(s, "\n") }
    write( paste(level, s, sep=" "), file=DEBUG_LOG$FILE, append=TRUE )
  }
  
  parseList <- function(l, br=3, indent="  ") {
    s <- indent
    counter <- 0
    n <- names(l)
    for (k in n) {
      v <- l[[k]]
      if (counter < br) {
        s <- paste(s, sprintf("%s=%s, ", k, v), sep="")
        counter <- counter + 1
      } else {
        s <- paste(s, "\n", indent, sprintf("%s=%s, ", k, v), sep="")
        counter <- 0
      }
    }
    paste(s, "\n", sep="")
  }
  
  list(FILE=FILE, LEVEL_WARN=LEVEL_WARN, LEVEL_INFO=LEVEL_INFO,
       LEVEL_DEBUG=LEVEL_DEBUG, LEVEL_TRACE=LEVEL_TRACE, LEVEL=LEVEL, CONSOLE=CONSOLE,
       INDENT=INDENT, OUTPUT=OUTPUT, parseList=parseList)
})()


TRACE <- function(..., i=0) {
  if (DEBUG_LOG$LEVEL >= DEBUG_LOG$LEVEL_TRACE) {
    DEBUG_LOG$OUTPUT("TRACE:", DEBUG_LOG$INDENT[i+1], sprintf(...))
  }
}

DEBUG <- function(..., i=0) {
  if (DEBUG_LOG$LEVEL >= DEBUG_LOG$LEVEL_DEBUG) {
    DEBUG_LOG$OUTPUT("DEBUG:", DEBUG_LOG$INDENT[i+1], sprintf(...))
  }
}

INFO <- function(..., i=0) {
  if (DEBUG_LOG$LEVEL >= DEBUG_LOG$LEVEL_INFO) {
    DEBUG_LOG$OUTPUT("INFO :", DEBUG_LOG$INDENT[i+1], sprintf(...))
  }
}

WARN <- function(..., i=0) {
  if (DEBUG_LOG$LEVEL >= DEBUG_LOG$LEVEL_WARN) {
    DEBUG_LOG$OUTPUT("WARN :", DEBUG_LOG$INDENT[i+1], sprintf(...))
    warning(sprintf(...))
  }
}

ERROR <- function(..., i=0) {
  if (DEBUG_LOG$LEVEL >= DEBUG_LOG$LEVEL_ERROR) {
    DEBUG_LOG$OUTPUT("ERROR:", DEBUG_LOG$INDENT[i+1], sprintf(...))
    warning(sprintf(...))
  }
}

FATAL <- function(..., i=0) {
  if (DEBUG_LOG$LEVEL >= DEBUG_LOG$LEVEL_FATAL) {
    DEBUG_LOG$OUTPUT("FATAL:", DEBUG_LOG$INDENT[i+1], sprintf(...))
    warning(sprintf(...))
  }
}

CAT <- function(..., i=0) {
  if (DEBUG_LOG$LEVEL >= DEBUG_LOG$LEVEL_INFO) {
    DEBUG_LOG$OUTPUT("INFO :", DEBUG_LOG$INDENT[i+1], sprintf(...), console=FALSE)
  }
  cat(paste(DEBUG_LOG$INDENT[i+1], sprintf(...), sep=""), "\n")
}


# _____________________________________________________________________________

# ** UP (Utility functions for plots) ** ----

UP <- (function(){
  
  col32_ <- function(col, alpha) {
    rgb_ = col2rgb(col, alpha=TRUE)
    ifelse(alpha >= 0,
           rgb(rgb_[1], rgb_[2], rgb_[3], alpha, maxColorValue=255),
           rgb(rgb_[1], rgb_[2], rgb_[3], maxColorValue=255))
  }
  
  list(
    # Add alpha level to a color.
    #   col - name or hex (e.g. red, #123456)
    #   alpha - number from 0 to 255 [optional]
    col32 = function(cols, alphas) {
      if (missing(alphas)) { alphas = -1 }
      mapply(col32_, cols, alphas)
    },
    
    errorBar = function(x0, y0, x1, y1, col, length) {
      arrows(x0, y0, x1, y1, angle=90, code=3, length=length, col=col)
    }, 
    errorBarX = function(x, y, err, col, length=0.1) {
      UP$errorBar(x-err, y, x+err, y, col, length)
    },
    errorBarY = function(x, y, err, col, length=0.1) {
      UP$errorBar(x, y-err, x, y+err, col, length)
    }
  )
})()





# ****************************************
# ** FCP (Flow in a circular pipe) ** ----
# ****************************************

FCP <- (function(){
  
  fD.colebrook_ <- function(roughness, D, Re, tol=1e-8, itMax=10, warn=TRUE) {
    core_ <- function(roughness, D, Re) {
      if (Re <= 4000 && (warn == TRUE)) { warning("Re <= 4000 !") }
      
      fun <- function(fD) {
        (1 / sqrt(fD)) + 2 * log10( roughness / D / 3.71 + 2.51 / Re / sqrt(fD))
      }
      
      # Derivative of fun().
      dFun <- function(fD) {
        - fD^(-3/2) * (1/2 + 2.51 / log(10) / (2.51 / Re / sqrt(fD) + roughness / 3.71 / D) / Re)
      }
      
      # Newton-Raphson method
      fD_0 <- FCP$fD.Blasius(Re)
      fD_n <- UTIL$newtonRaphson(fun, dFun, fD_0, tol=tol, itMax=itMax) 
      fD_n
    }
    mapply(core_, roughness, D, Re)
  }
  
  list(
    
    # Darcy-Weisbach equation (for the calculation of dP per unit length)
    DarcyWeisbach = function(fD, density, velocity, D) {
      fD * density * (velocity ^ 2) / (2 * D)
    },
    
    WallShearStress = function(D, dp_dL) {
      D / 4 * dp_dL
    },
    
    frictionalVelocity = function(WallShearStress, density){
      sqrt(WallShearStress / density)
    },
    
    # Law of the wall
    #   y+ < 5:         Viscous sublayer
    #   5 < y+ < 30~70: Buffer layer
    yPlus = function(frictionalVelocity, y, kinematicViscosity) {
      frictionalVelocity * y / kinematicViscosity    # y: distance from the wall 
    },
    
    #uPlus <- function(velocity, frictionalVelocity) {
    #  velocity / frictionalVelocity
    #},
    
    # Darcy friction factor
    fD.laminar = function(Re) { 64 / Re },
    fD.Blasius = function(Re, C=0.3164) { C / (Re ^ 0.25) },
    fD.Colebrook = fD.colebrook_,
    
    # Fanning friction factor
    ff.laminar = function(Re) { 16 / Re },
    ff.Blasius = function(Re, C=0.0791) { C / (Re ^ 0.25) },
    ff.Colebrook = function(roughness, D, Re, tol=1e-8, warn=TRUE){
      fD.colebrook_(roughness, D, Re, tol=tol, warn=warn) / 4
    },
    
    
    # Borda-Carnot's formula (Energy loss)
    #   dE = lossCoefficient * density * (vout - vin)^2 / 2
    #     vin: Flow velocity before expansion
    #     Ain: Cross sectional area before expansion
    #     Aout: Cross sectional area after expansion
    #     density: Fluid density
    BordaCarnot = function(vin, Ain, Aout, density){
      (1 - (Ain / Aout))^2 * density * vin^2 / 2
    }, 
    
    # Borda-Carnot (Head loss)
    #   vin : Flow velocity before expansion
    #   Ain: cross sectional area before expansion
    #   Aout: cross sectional area after expansion
    BordaCarnotHead = function(vin, Ain, Aout){
    	(1 - (Ain / Aout))^2 * (vin^2) / (2 * g)
    },
    
    
    # Drift Flux Model.
    driftFluxModel = function(C0, Vd, vsg, vsl) {
      vmix <- vsg + vsl
      vg <- C0 * vmix + Vd
      voidRatio <- vsg / vg
      
      Hl <- 1 - voidRatio
      vl <- vsg / Hl
      
      c(vg=vg, vl=vl, Hl=Hl)
    }
    
  )
  
})()



# Heat Exchanger 


# ******************************************************************
# ** LMTD (logarithmic mean temperature difference) ** ----
#
#   - Positive values mean heat flow from line-1 to line-2 (T1 > T2) 
#   - Nagative values mean heat flow from line-2 to line-1 (T2 < T1)
# ******************************************************************

LMTD = (function() {
  
  assert = function(expr, msg, T1_in, T1_out, T2_in, T2_out) {
    if (expr == FALSE) {
      errorMsg = sprintf("%s T1_out=%.2f, T2_out=%.2f, T1_in=%.2f, T2_in=%.2f", msg, T1_in, T1_out, T2_in, T2_out)
      stop(errorMsg)
    }
  }

  list(

    # Meaning of the symbols:
  	#   dT   in         out
  	#
  	#   T1   in ------> out
  	#   T2  out <------ in
  	counterCurrent = function(T1_in, T1_out, T2_in, T2_out) {
  		# Error check of inputs
  		if (T1_in > T2_out) {
  			assert(T1_out > T2_in, "Temperature of line-1 must be higher at both ends (T1 > T2).", T1_in, T1_out, T2_in, T2_out)
  			assert(T1_in > T1_out, "Fluid of line-1 must be cooled down (T1_in > T1_out) because of 'T1_in > T2_out'.", T1_in, T1_out, T2_in, T2_out)
  			assert(T2_in < T2_out, "Fluid of line-2 must be warmed up (T2_in < T2_out) because of 'T1_in > T2_out'.", T1_in, T1_out, T2_in, T2_out)
  		} else if (T1_in < T2_out) {
  			assert(T1_out < T2_in, "Temperature of line-1 must be lower at both ends (T1 < T2).", T1_in, T1_out, T2_in, T2_out)
  			assert(T1_in < T1_out, "Fluid of line-1 must be warmed up (T1_in < T1_out) because of 'T1_in < T2_out'.", T1_in, T1_out, T2_in, T2_out)
        assert(T2_in > T2_out, "Fluid of line-2 must be cooled down (T2_in > T2_out) because of 'T1_in < T2_out'.", T1_in, T1_out, T2_in, T2_out)
      } else {
        # some check may be required.
      }
    
      dT_in  <- T1_in - T2_out
      dT_out <- T1_out - T2_in
      (dT_in - dT_out) / log(dT_in / dT_out)
    },
  
    # Meaning of the symbols:
    #   dT   in         out
    #
    #   T1   in ------> out
    #   T2   in ------> out
    coCurrent = function(T1_in, T1_out, T2_in, T2_out) {
      # Error check of inputs
      if (T1_in > T2_in) {
        assert(T1_out > T2_out, "Temperature of line-1 must be higher at both ends (T1 > T2).", T1_in, T1_out, T2_in, T2_out)
        assert(T1_in > T1_out, "Fluid of line-1 must be cooled down (T1_in > T1_out) because of 'T1_in > T2_in'.", T1_in, T1_out, T2_in, T2_out)
        assert(T2_in < T2_out, "Fluid of line-2 must be warmed up (T2_in < T2_out) because of 'T1_in > T2_in'.", T1_in, T1_out, T2_in, T2_out)
      } else if (T1_in < T2_in) {
        assert(T1_out < T2_out, "Temperature of line-1 must be lower at both ends (T1 < T2).", T1_in, T1_out, T2_in, T2_out)
        assert(T1_in < T1_out, "Fluid of line-1 must be warmed up (T1_in < T1_out) because of 'T1_in < T2_out'.", T1_in, T1_out, T2_in, T2_out)
        assert(T2_in > T2_out, "Fluid of line-2 must be cooled down (T2_in > T2_out) because of 'T1_in < T2_out'.", T1_in, T1_out, T2_in, T2_out)
      } else {
        assert(T1_out == T2_out, "Temperature of both fluids must be same at outlet (T1_out == T2_out).", T1_in, T1_out, T2_in, T2_out)
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

})()




## ** IdealGas ** ----

IdealGas <- (function() {

  Smix   <- function(n, x_i) { - n * R * sum(x_i * log(x_i)) }

  list(
    # EOS: PV = nRT
    P = function(V, n, T) { n * R * T / V },    # Pa
    V = function(P, n, T) { n * R * T / P },    # m3
    n = function(P, V, T) { P * V / R / T },    # mol
    T = function(P, V, n) { P * V / n / R },    # K
    
    # Raoult's law
    Raoult.Pi = function(x_i, Psat_i) { x_i * Psat_i },
    Raoult.BP = function(x_i, Psat_i) { sum( x_i * Psat_i ) },     # Bubble Point
    Raoult.DP = function(x_i, Psat_i) { 1 / sum(x_i * Psat_i) },   # Dew Point
    
    # Mixing Entropy and Mixing Gibbs Energy
    Smix = Smix,
    Gmix = function(n, x_i, T) { - n * T * Smix(n, x_i) }

  )
  
})()




m <- list(
  
  # mode()
  
  
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Greatest Common Divisor & Least Common Multiple ----
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  gcd = (function(){
    gcd_ <- function(a, b){
      while(a %% b != 0){
        tmp <- b
        b <- a %% b
        a <- tmp
      }
      return(b)
    }
    function(...){ Reduce(gcd_, c(...)) }
  })(),
  
  lcm = (function(){
    lcm_ <- function(a, b){
      a * b / MATH$gcd(a,b)
    }
    function(...){ Reduce(lcm_, c(...)) }
  })()
  
)

