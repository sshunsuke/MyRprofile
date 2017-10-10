




hydraulicDiameter <- function(A, S) {
  4 * A / S
}

cpf.DarcyWeisbach <- function(fD, density, velocity, D) {
  fd * density * (velocity ^ 2) / (2 * D)
}

# =============================================================================
# Friction Factor.
# =============================================================================

# Darcy Friction Factor
cpf.fd <- list(
  Colebrook = function(roughness, D, Re, interval=c(0, 5)) {
    cpf.ff$Colebrook(roughness, D, Re, interval=interval)
  },
  Blasius = function(Re, C=0.3164) { C / (Re ^ 0.25) },
  laminar = function(Re) { 64 / Re }
)

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
    
    Colebrook = colebrook_,
    Blasius = function(Re, C=0.0791) { C / (Re ^ 0.25) },
    laminar = function(Re) { 16 / Re }
  )
})()


# =============================================================================
# Drift Flux Model.
# =============================================================================

cpf.dfm <- (function(){
  
  
  
})()

