# =============================================================================
# Constant values
# =============================================================================

g <- 9.8                      # Gravity Acceleration (m/s2)
R <- 8.3144621                # Gas Constant (J/K-mol)
# R <- 8.314471               # Moldover et al. (1988)
# Na <- 6.022140857 * 10^23     # Avogadro constant
kB <- 1.38064852 * 10^(-23)   # Boltzmann constant (J/K)

# STP (Standard Temperature and Pressure) - old 
Pstp <- 100000    # (Pa)
Tstp <- 273.15    # (K)

# SATP (Standard Ambient Temperature and Pressure)
Psatp <- 101325    # (Pa)
Tsatp <- 298.15    # (K)

# NTP (Normal Temperature and Pressure)
Pntp <- 101325    # (Pa)
Tntp <- 293.15    # (K)


# =============================================================================
# Unit Conversion
# =============================================================================
if (exists('inch2m') == TRUE) { detach(UC) }

UC <- list(
  rad2deg = function(rad) { rad * 180 / pi},
  deg2rad = function(deg) { deg * pi / 180},
  
  inch2m = function(inch) { inch * 0.0254 },
  m2inch = function(m) { m / 0.0254 },
  ft2m   = function(ft) { ft * 0.3048 },
  m2ft   = function(m) { m / 0.3048 },
  
  K2C = function(K) { K - 273.15 },
  C2K = function(C) { C + 273.15 },
  F2C = function(F) { (F - 32) * 5 / 9 },
  C2F = function(C) { (C * 9/5) + 32 },
  F2R = function(F) { F + 459.67 },
  R2F = function(R) { R - 459.67 },
  
  psi2Pa = function(psi) { psi * 6894.76 }
)

if (exists('inch2m') == FALSE) { attach(UC) }


# =============================================================================
# Properties of simple shapes.
# =============================================================================
circle <- list(
  area = function(r) { r * r * pi },
  circumference = function(r) { 2 * r * pi }
)

sphere <- list(
  volume      = function(r) { (4 * pi * r^3) / 3 },
  surfaceArea = function(r) { 4 * pi * r^2 }
)

elipsoid <- list(
  volume      = function(a,b,c) { 4/3 * pi * a * b * c }
)

# Sauter Mean Diameter (vD: Vector of diameters)
SMD <- function(vD) { sum(vD^3) / sum(vD^2) }


# =============================================================================
# Dimensionless Number
# =============================================================================
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

    
    
    
    # -------------------------------------------------------------------------
    # Ratios between diffusions of momentum, heat or mass
    # -------------------------------------------------------------------------
    
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
      k / (density * specificHeat * diffusionCoefficient)
    },
    
    
    # -------------------------------------------------------------------------
    # Ratio of different forces
    # -------------------------------------------------------------------------
    
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

# =============================================================================
# IO clipboard
# =============================================================================
rcbmat <- function(header, ...) {
  if (missing(header)) { header = FALSE }
  utils::read.table(file="clipboard", header=header, ...)
}

wcbmat <- function(data, header, sep="\t", row.names=FALSE, qmethod="double",
                   col.names=ifelse(header && row.names, NA, header), ...) {
  if (missing(header)) { header = FALSE }
  utils::write.table(data, file="clipboard-128", sep=sep, row.names=row.names,
                     col.names=col.names, qmethod=qmethod, ...)
}

# =============================================================================
# Vector
# =============================================================================
v.indexClosestValue <- function(vec, x) {
	diff <- abs(vec-x)
	which( diff == min(diff) )
}

v.closestValue <- function(vec, x) {
	vec[ v.indexClosestValue(vec, x) ]
}



# =============================================================================
# Matrix
# =============================================================================
m.crev <- function(mat) { mat[,ncol(mat):1] }
m.rrev <- function(mat) { mat[nrow(mat):1,] }

# =============================================================================
# Data Frame
# =============================================================================
df.orderBy <- function(df, colname, decreasing=FALSE) {
  if (missing(colname)) { stop("'colname' is not specified.") }
  df[order(df[,colname], decreasing=decreasing),]
}

# =============================================================================
# Functions for debug log.
# =============================================================================
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



# =============================================================================
# Flow in a circular pipe.
# =============================================================================

FCP <- (function(){
  
  fD.colebrook_ <- function(roughness, D, Re, tol=1e-8, warn=TRUE) {
    core_ <- function(roughness, D, Re) {
      if (Re <= 4000 && (warn == TRUE)) { warning("Re <= 4000 !") }
      
      fun <- function(fD) {
        (1 / sqrt(fD)) + 2 * log10( roughness / D / 3.71 + 2.51 / Re / sqrt(fD))
      }
      
      # Derivative of fun().
      dFun <- function(fD) {
        - fD^(-3/2) * (1/2 + 2.51 / log(10) / (2.51 / Re / sqrt(fD) + roughness / 3.71 / D) / Re)
      }
      
      fD_n <- FCP$fD.Blasius(Re)
      d <- fun(fD_n) / dFun(fD_n)
      
      # Newton-Raphson method
      while (abs(d) >= tol) {
        d <- fun(fD_n) / dFun(fD_n)
        fD_n <- fD_n - d
      }
      
      fD_n
    }
    mapply(core_, roughness, D, Re)
  }
  
  list(
    DarcyWeisbach = function(fD, density, velocity, D) {
      fd * density * (velocity ^ 2) / (2 * D)
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
    }
  )
  
})()


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



