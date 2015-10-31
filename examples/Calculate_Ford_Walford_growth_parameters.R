
###################
# NEW PROPOSAL
#1.  Find VGBF "k" parameter from fishbase or elsewhere (I'll call it k_L for growth coefficient K in the length VBGF)
#2.  Find natural mortality rate "M".  If unavailable, calculate from simplest life history, M = k * 1.6
#3.  Find Linf from fishbase or whatever
#4.  Find L0 (or t0) from for VBGF in length.  if unavailable, assume L0=t0=0
#5.  Find weight at age parameters a and b, such that W(a) = a*L(a)^b
#6.  Find age at recruitment "r", and weight at this age, "w_r"
#6b. If age at recruitment is missing, use age-at-maturity, rounded to nearest age (or up to 1 if <1)
#6c. If age at maturity is missing, calculate likely age at maturity from life history theory, r = log( (3*k+M)/M ) / K, rounded to nearest age (or up to 1 if <1)
#7. Determine average weight in your data set.  I'll call this W(abar)
#
# Derived values
#A. Calculate asymptotic maximum weight Winf = a*Linf^b
#B. Calculate age corresponding to W(abar) from VBGF in length, given k, Linf, L0, a, and b, rounded to nearest age (or up to 1 if <1).  I'll call this abar
#C. Calculate weight for abar+1, W(abar+1) from VBGF in length, given k, Linf, L0, a, and b 
#D. Calculate rho such that annual growth in Ford-Walford form matches the difference between W(abar+1) and W(abar).  rho = (Winf - W(abar+1))/(Winf - W(abar))
#E. Calculate a_g = Winf * (1-rho)
#F. Calculate W_r, weight at r, 
#####################################################


Ford_Walford_Params = function( k, M=1.6*k, Linf, L0=0, a, b=3, r=round(log((3*k+M)/M)/k), Wbar=a*Linf^b/2 ){
  # maximum weight
  Winf = a*Linf^b
  amax = log(0.01) / -M
  # VBGF
  a_set = seq(0,amax,by=0.1)
  L_a = Linf - (Linf-L0)*exp(-k*a_set)
  abar = round(a_set[which.min( abs(a*L_a^b - Wbar) )])
  W_abar = a * (Linf - (Linf-L0)*exp(-k*abar))^b
  W_abar_plusone = a * (Linf - (Linf-L0)*exp(-k*(abar+1)))^b  
  rho = (Winf - W_abar_plusone)/(Winf - W_abar)
  a_g = Winf * (1-rho)
  # Calculate Ford-Walford growth curve
  a_set = 1:amax
  W_a = rep(NA,length(a_set))
  W_a[abar] = W_abar
  for(age in (abar+1):amax) W_a[age] = rho*W_a[age-1] + a_g  
  for(age in (abar-1):1) W_a[age] = (W_a[age+1] - a_g) / rho
  # 
  W_r = W_a[r]
  # Plot results
  plot( x=0:amax, y=a*(Linf - (Linf-L0)*exp(-k*0:amax))^b, type="l", lwd=2, ylab="Weight (g)", xlab="Age (years)" )
  lines( x=r:amax, W_a[r:amax], col="red", lwd=2)
  abline( v=c(r,abar), lty=c("solid","dashed"))
  # 
  if( W_a[r]<0 ) stop("Ford-Walford weight at age for age at recruitment is less than zero, please recompute")
  # Return stuff
  Return = list("Winf"=Winf, "rho"=rho, "a_g"=a_g, "W_r"=W_r)
  return(Return)
}

Ford_Walford_Params( k=0.01, Linf=30, a=0.01 )

