# =============================================================================
# Constant values
# =============================================================================
# Gravity Acceleration (m/s2)
g <- 9.8

# Gas Constant (J/K-mol)
R <- 8.3144621
# R <- 8.314471    # Moldover et al. (1988)

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
  inch2m = function(inch) { inch * 0.0254 },
  m2inch = function(m) { m / 0.0254 },
  ft2m   = function(ft) { ft * 0.3048 },
  m2ft   = function(m) { m / 0.3048 },
  
  K2C = function(K) { K - 273.15 },
  C2K = function(C) { C + 273.15 },
  F2C = function(F) { (F - 32) * 5 / 9 },
  
  psi2Pa = function(psi) { psi * 6894.76 }
)

if (exists('inch2m') == FALSE) { attach(UC) }


# =============================================================================
# Dimensionless Number
# =============================================================================
if (exists('Reynolds') == TRUE) { detach(DN) }

DN <- list(
  
  Weber = function(density, v, L, surfaceTension) {
    density * (v^2) * L / surfaceTension
  },
  
  # specificHeat: J/Kg-K, viscosity: N-s/m2, thermal conductivity: W/m-K
  Prandtl = function(specificHeat, viscosity, thermalConductivity) {
    (specificHeat * viscosity) / thermalConductivity
  },
  
  Reynolds = function(density, v, L, viscosity) {
    density * v * L / viscosity
  }
  
)

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
# Matrix
# =============================================================================
m.crev <- function(mat) { mat[,ncol(mat):1] }
m.rrev <- function(mat) { mat[nrow(mat):1,] }

# =============================================================================
# Data Frame
# =============================================================================
df.orderBy <- function(df, colname) { df[order(df[,colname]),] }

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





circle <- list(
  area = function(r) { r * r * pi },
  circumference = function(r) { 2 * r * pi }
)

sphere <- list(
  volume      = function(r) { (4 * pi * r^3) / 3 },
  surfaceArea = function(r) { 4 * pi * r^2 }
)





