




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
# =============================================================================

cpf.dfm <- function(C0, Vgj, vsg, vsl) {
  vmix <- vsg + vsl
  vg <- vmix + Vgj
  Hl <- 1 - (vsg / vg)
  vl <- vsg / Hl
  
  c(Hl = Hl, vg=vg, vl=vl)
}

cpf.dfm.C0 <- list(
  bubble = function(densityG, densityL) {
    
  }
  
  
)

