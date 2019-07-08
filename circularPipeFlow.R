




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

