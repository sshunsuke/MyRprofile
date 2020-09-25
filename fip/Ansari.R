
Ansari <- list()
Ansari$Util <- list()


# Slug flow ----

Ansari$Slug <- function(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, tol=1e-8) {
	vmix = vsG + vsL
	
	# Bubble Rise Velocity - 4.160
	vbr = 1.53 * ((g * surfaceTension * (densityL - densityG)) / densityL^2) ^ 0.25
	
	# Void ration (and holdup) of liquid slug - 4.194
	Hgls = vsG / (0.425 + 2.65 * vmix)
	Hlls = 1 - Hgls
	
	# Taylor-bubble-rise Velocity - 4.190
	vtb = 1.2 * vmix + 0.35 * ((g * D * (densityL - densityG)) / densityL^2) ^ 0.5
	
	# Gas velocity in liquid slug - 4.191
	vgls = 1.2 * vmix + vbr * Hlls^0.5
	
	# Holdup (and void ratio) in Taylor bubble section (4.196-199)
	Hltb = Ansari$SlugHltb(D, vmix, Hgls, vtb, vgls, tol=tol)
	Hgtb = 1 - Hltb
	
	# Velocity of falling liquid film in Taylor-bubble section (4.193)
	vltb = - 9.916 * ( g * D * (1 - sqrt(Hgtb)) )^0.5
	
	# Mass balance of liquid with respect to Taylor-bubble section (4.188)
	vlls = vtb - (vtb - vltb) * (Hltb / Hlls)
	#Vlls = (vmix - Vgls * (1 - Hlls)) / Hlls
	
	# Mass balance of gas with respect to Taylor-bubble section (4.189)
	vgtb = vtb - (vtb - vgls) * (Hgls / Hgtb)
	
	# --- --- --- --- --- --- --- --- 
	
	# Proportion of Slug
	filmThickness = (vltb^2) / (196.7 * g)
	
	# Ratio between Taylor bubble and liquid slug (4.185-4.187)
	Ltb_Lu <- (vsL - vlls * Hlls) / (vltb * Hltb - vlls * Hlls)
	Lls_Lu <- 1 - Ltb_Lu
	
	Hlu <- Ltb_Lu * Hltb + Lls_Lu * Hlls
	
	data.frame(Hgls, Hlls, vgls, vlls, Hgtb, Hltb, vgtb, vltb, Hlu, filmThickness)
}


Ansari$SlugHltb <- function(D, vmix, Hgls, vtb, vgls, tol=1e-8) {
	# 4.196, 4.197
	fun <- function(Hltb) {
		A = Hgls * (vtb - vgls) + vmix 
		(9.916 * sqrt(g * D) * (1 - sqrt(1- Hltb))^0.5 * Hltb) - vtb * (1 - Hltb) + A
	}
	
	# 4.198
	dFun <- function(Hltb) {
		a <- 1 - sqrt(1- Hltb)
		vtb + (9.916 * sqrt(g * D) * a^0.5) + Hltb / (4 * sqrt(1 - Hltb * a))
	}
	
	# 4.199
	Hltb_j <- 0.15
	d <- fun(Hltb_j) / dFun(Hltb_j)
	
	while (abs(d) >= tol) {
		d <- fun(Hltb_j) / dFun(Hltb_j)
		Hltb_j <- Hltb_j - d
	}
	
	return(Hltb_j)
}



# Util ----

# Bubble Rise Velocity.
# 
# References:
#   Harmathy, T. Z. (1960).
#   Velocity of large drops and bubbles in media of infinite or restricted extent.
#   AIChE Journal, 6(2), 281-288. doi:10.1002/aic.690060222
Ansari$Util$VbrHarmathy <- function(densityG, densityL, surfaceTension) {
	1.53 * ((g * surfaceTension * (densityL - densityG)) / densityL^2) ^ 0.25
}



# ************* ----

# test ----

# flowProperties(P, T, diameter, Qgstd, Ql, ipDiameter=0, inclination=pi/2)

#fp_ <- flowProperties(200*1000, 273.15+10, inch2m(5.9), 24000/24/3600, 200/24/3600, inclination=pi/2)
#dPcalc(fp_, 10, FR$SLUG)

# * thickness: 0.005244 (m),  Lls_Lu: 0.126009,  Ltb_Lu: 0.873991 
# * Hlls: 0.635680,  Hltb: 0.135073,  Hlu: 0.198155,   
# - Vmix: 8.279308 (m/sec),  Vtb: 10.184623 (m/sec) 
# - Vgls: 10.134058 (m/sec),  Vgtb: 10.163324 (m/sec),  Vlls: 7.216317 (m/sec),  Vltb: -3.784761 (m/sec) 
# - densityU: 205 (kg/m3),  densityLS: 654.0 (kg/m3),  viscosityLS: 0.000894 (Pa-s) 
# - ReLS: 907410.0,  fLS: 0.003999,  SSls: 89.6 (N/m2),  SSwf: -50.3 (N/m2) 
# * dP_H: 2007.0 (Pa/m),  dP_F: -872.9 (Pa/m),  dP_A: 0.0 (Pa/m) 
# * diffP: 11340.7 (Pa),  dZ: 10.0 (m) 

#Ansari$Slug <- function(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, tol=1e-8) 


Ansari$Slug(8.148071, 0.1312366, 0.14986, 1.368176,  1027.98, 1.065058e-05, 0.001400565, 0.07422209)




