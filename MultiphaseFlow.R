

Ansari <- list()
Ansari$Util <- list()


Ansari$Slug <- function(vsG, vsL, D, densityG, densityL, viscosityG, viscosityL, surfaceTension, tol=1e-8) {
	Vmix = vsG + vsL
	
	# Bubble Rise Velocity - 4.160
	Vbr = 1.53 * ((g * surfaceTension * (densityL - densityG)) / densityL^2) ^ 0.25
	
	# Void ration (and holdup) of liquid slug - 4.194
	Hgls = vsG / (0.425 + 2.65 * Vmix)
	Hlls = 1 - Hgls
	
	# Taylor-bubble-rise Velocity - 4.190
	Vtb = 1.2 * Vmix + 0.35 * ((g * D * (densityL - densityG)) / densityL^2) ^ 0.5
	
	# Gas velocity in liquid slug - 4.191
	Vgls = 1.2 * Vmix + Vbr * Hlls^0.5
	
	# Holdup (and void ratio) in Taylor bubble section (4.196-199)
	Hltb = Ansari$SlugHltb(D, vmix, Hgls, vtb, vgls, tol=tol)
	Hgtb = 1 - Hltb
	
	# Velocity of falling liquid film in Taylor-bubble section (4.193)
	Vltb = - 9.916 * ( g * D * (1 - sqrt(Hgtb)) )^0.5
	
	# Mass balance of liquid with respect to Taylor-bubble section (4.188)
	Vlls = Vtb - (vtb - Vltb) * (Hltb / Hlls)
	#Vlls = (Vmix - Vgls * (1 - Hlls)) / Hlls

	# Mass balance of gas with respect to Taylor-bubble section (4.189)
	vgtb = vtb - (vtb - Vgls) * (Hgls / Hgtb)
	
	# ----------------------------------------
	
	# Proportion of Slug
	filmThickness = (vltb^2) / (196.7 * g)
	
	# Ratio between Taylor bubble and liquid slug (4.185-4.187)
	Ltb_Lu <- (vsL - Vlls * Hlls) / (Vltb * Hltb - Vlls * Hlls)
	Lls_Lu <- 1 - Ltb_Lu
	
	Hlu <- Ltb_Lu * Hltb + Lls_Lu * Hlls
	
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








# Bubble Rise Velocity.
# 
# References:
#   Harmathy, T. Z. (1960).
#   Velocity of large drops and bubbles in media of infinite or restricted extent.
#   AIChE Journal, 6(2), 281-288. doi:10.1002/aic.690060222
Ansari$Util$VbrHarmathy <- function(densityG, densityL, surfaceTension) {
	1.53 * ((g * surfaceTension * (densityL - densityG)) / densityL^2) ^ 0.25
}








dPcalc.slug.V.dP <- function(fp, dZ) {
	INFO("START:  dPCalc.slug.V.dP()")
	
	dDensity <- fp$densityL - fp$densityG
	
	# Holdup in liquid slug - Hasan & Kabir (1988)
	#Vtb <- dPcalc.correlation.slug.Vtb.HasanKabir(fp, dDensity)
	#Vtb <- dPcalc.correlation.slug.Vtb.Bendiksen(fp)
	Vtb <- 1.2 * fp$Vmix + dPcalc.correlation.Vbr.Harmathy(fp, dDensity)
	
	# Holdup in liquid slug - Sylvester (1987) or Barnea & Brauner (1985)
	Hlls <- dPcalc.correlation.slugV.Hlls(fp)
	
	# Gas velocity in liquid slug - Ansari (1994)
	Vgls <- dPcalc.correlation.slug.Vgls.Ansari(fp, dDensity, Hlls)
	#Vgls <- dPcalc.correlation.slug.Vgls.Barnea(fp, dDensity)
	
	Hlu_ <- (Vtb * Hlls + Vgls * (1 - Hlls) - fp$Vsg) / Vtb
	
	DEBUG("  Hlls: %f,  Hlu: %f,  Vgls: %f (m/sec), Vtb: %f (m/sec)", Hlls, Hlu_, Vgls, Vtb)
	DEBUG("  - - - - -")
	
	fun <<- function(thickness_) {
		DEBUG("  fun():")
		DEBUG("    thickness_ = %f (m)", thickness_)
		tbp_ <- dPcalc.slug.V.properties(fp, thickness_, Hlls, Vgls, Vtb)
		if (DEBUG_LOG$LEVEL >= DEBUG_LOG$LEVEL_TRACE) {
			TRACE( dPcalc.slug.V.flowProperties.toString(tbp_) )
		}
		
		#
		if (tbp_$Vgtb < Vgls) {
			DEBUG("    tbp_$Vgtb < Vgls  [tbp_$Vgtb=%f, Vgls=%f]", tbp_$Vgtb, Vgls)
			return (-Inf)
		} 
		a <- tbp_$SSwf * tbp_$Sl / tbp_$Al
		b <- tbp_$SSi * tbp_$Si * (1/tbp_$Al + 1/tbp_$Ag)
		c <- (fp$densityL - fp$densityG) * g * sin(fp$inclination)
		DEBUG("    a=%f, b=%f, c=%f, ret=%f", a, b, c, a-b+c)
		a - b + c
	}
	
	thicknessList <- seq(dPcalc.configuration$slugV.minThickness, fp$De/2-0.001, length.out=50)
	retFun <- sapply(thicknessList, fun)
	maxFun <- max(retFun)
	
	# 
	if (maxFun < 0) {
		thicknessList <- seq(thicknessList[which(retFun == max(retFun))], thicknessList[which(retFun == max(retFun))+1], length.out=50)
		retFun <- sapply(thicknessList, fun)
		#cat(retFun)
		maxFun <- max(retFun)
	}
	
	range <- c(dPcalc.configuration$slugV.minThickness, thicknessList[which(retFun == max(retFun))])
	
	TRACE("  - - -")
	TRACE("  maxFun=%f, range[2]=%f (m)", maxFun, range[2])
	TRACE("* Start finding solution...")
	TRACE("*   Thickness Range: %f ~ %f (m)", range[1], range[2])
	
	f <- uniroot(fun, range, tol=1e-10)
	
	DEBUG("  - - - - -")
	
	thickness <- f$root
	tbp <- dPcalc.slug.V.properties(fp, thickness, Hlls, Vgls, Vtb)
	
	Ltb_Lu <- (fp$Vsl - tbp$Vlls * Hlls) / (tbp$Vltb * tbp$Hltb - tbp$Vlls * Hlls)
	Lls_Lu <- 1 - Ltb_Lu
	
	Hlu <- Ltb_Lu * tbp$Hltb + Lls_Lu * Hlls
	
	#cat( sprintf("%f, Hlu_: %f\n", Hlu, Hlu_))
	#cat( sprintf("%f \n", Ltb_Lu))
	
	#Lls_Lu <- (Hlu - tbp$Hltb) / (Hlls - tbp$Hltb)
	#Ltb_Lu <- 1 - Lls_Lu
	
	INFO("* thickness: %f (m),  Lls_Lu: %f,  Ltb_Lu: %f", thickness, Lls_Lu, Ltb_Lu)
	INFO("* Hlls: %f,  Hltb: %f,  Hlu: %f,  ", Hlls, tbp$Hltb, Hlu)
	INFO("- Vmix: %f (m/sec),  Vtb: %f (m/sec)", fp$Vmix, Vtb, Vgls)
	INFO("- Vgls: %f (m/sec),  Vgtb: %f (m/sec),  Vlls: %f (m/sec),  Vltb: %f (m/sec)", Vgls, tbp$Vgtb, tbp$Vlls, tbp$Vltb)
	
	densityU <- Hlu * fp$densityL + (1 - Hlu) * fp$densityG
	densityLS <- Hlls * fp$densityL + (1 - Hlls) * fp$densityG
	viscosityLS <- Hlls * fp$viscosityL + (1 - Hlls) * fp$viscosityG
	
	INFO("- densityU: %.f (kg/m3),  densityLS: %.1f (kg/m3),  viscosityLS: %f (Pa-s)", densityU, densityLS, viscosityLS)
	
	ReLS <- densityLS * fp$Vmix * fp$De / viscosityLS
	fLS <- dPcalc.fFrictionFactor(ReLS, fp$De)
	SSls <- fLS * densityLS * fp$Vmix^2 / 2
	
	INFO("- ReLS: %.1f,  fLS: %f,  SSls: %.1f (N/m2),  SSwf: %.1f (N/m2)", ReLS, fLS, SSls, tbp$SSwf)
	
	dP_Ftb <- tbp$SSwf * (fp$De * pi) / tbp$Ap * Ltb_Lu
	dP_Fls <- SSls * (fp$De * pi) / tbp$Ap * Lls_Lu
	
	dP_H <- densityU * g * sin(fp$inclination)
	dP_F <- dP_Ftb + dP_Fls
	dP_A <- 0   # (tbp$Vlls - tbp$Vltb) * (Vtb - tbp$Vlls) * fp$densityL * Hlu
	
	diffP <- (dP_H + dP_F + dP_A) * dZ
	
	INFO("* dP_H: %.1f (Pa/m),  dP_F: %.1f (Pa/m),  dP_A: %.1f (Pa/m)", dP_H, dP_F, dP_A)
	INFO("* diffP: %.1f (Pa),  dZ: %.1f (m)", diffP, dZ)
	INFO("END:  dPCalc.slug.V.dP()")
	
	#cat( (densityLS * g * sin(fp$inclination) + SSls * (fp$De * pi) / tbp$Ap) * Lls_Lu ); cat(" -- ")
	#cat( tbp$dP_lls * Ltb_Lu )
	
	list(diffP=diffP, Hl=Hlu,
			 dP_H=dP_H, dP_F=dP_F, dP_A=dP_A,
			 dP_Ftb=dP_Ftb, dP_Fls=dP_Fls, sfp=tbp, fun=fun)
}





dPcalc.slug.V.properties <- function(fp, thickness, Hlls, Vgls, Vtb) {
	Dg <- fp$De - thickness * 2    # Diameter of Taylor Bubble
	Sl <- fp$De * pi
	Si <- Dg * pi
	Ap <- (fp$De / 2)^2 * pi
	Ag <- (Dg / 2)^2 * pi
	Al <- Ap - Ag
	
	Hltb <- Al / Ap
	
	#TRACE("    thickness=%f(m), Dg=%f(m), Sl=%f(m), Si=%f(m)",thickness, Dg, Sl, Si)
	#TRACE("    Ap=%f(m2), Ag=%f(m2), Al=%f(m2), Hltb=%f", Ap, Ag, Al, Hltb)
	
	Vlls <- dPcalc.slug.Vlls(fp$Vmix, Hlls, Vgls)
	Vltb <- dPcalc.slug.Vltb(Hltb, Hlls, Vtb, Vlls)
	Vgtb <- dPcalc.slug.Vgtb(Hltb, fp$Vmix, Vltb)
	
	Dl <- 4 * Al / Sl 
	ReL <- (Dl * Vltb * fp$densityL) / fp$viscosityL  
	fL <- dPcalc.fFrictionFactor(abs(ReL), Dl)
	
	ReG <- (Dg * Vgtb * fp$densityG) / fp$viscosityG
	fG <- dPcalc.fFrictionFactor(abs(ReG), Dg, roughness=0)
	fI <- fG
	
	SSwf <- fL * fp$densityL * Vltb * abs(Vltb) / 2
	SSi  <- fI * fp$densityG * (Vgtb - Vltb) * abs(Vgtb - Vltb) / 2  #
	
	#TRACE("    Vlls=%f(m/s), Vltb=%f(m/s), Vgtb=%f(m/s)", Vlls, Vltb, Vgtb)
	#TRACE("    ReG=%.1f, ReL=%.1f, fG=%f, fL=%f", ReG, ReL, fG, fL)
	#TRACE("    SSwf=%f(N/m2), SSi=%f(N/m2)", SSwf, SSi)
	
	dP_lls <- fp$densityG * g * sin(fp$inclination) + SSi * Si / Ag
	dP_gls <- (SSwf * Sl - SSi * Si) / Al + fp$densityL * g * sin(fp$inclination)
	
	list(Dg=Dg, Ap=Ap, Ag=Ag, Al=Al, 
			 Hltb=Hltb, Vlls=Vlls, Vltb=Vltb, Vgtb=Vgtb,
			 ReL=abs(ReL), fL=fL, ReG=abs(ReG), fG=fG, fI=fI,
			 SSwf=SSwf, SSi=SSi, Sl=Sl, Si=Si,
			 dP_lls=dP_lls, dP_gls=dP_gls)
}


# Gomez (28)
dPcalc.slug.Vlls <- function(Vmix, Hlls, Vgls) {
	(Vmix - Vgls * (1 - Hlls)) / Hlls
}

# Gomez (27)
dPcalc.slug.Vltb <- function(Hltb, Hlls, Vtb, Vlls) {
	Vtb + (Vlls - Vtb) * Hlls / Hltb
}

# Gomez (29)
dPcalc.slug.Vgtb <- function(Hltb, Vmix, Vltb) {
	(Vmix - Vltb * Hltb) / (1 - Hltb)
}


