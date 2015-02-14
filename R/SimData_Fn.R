SimData_Fn <-
function(n_t, CV_C, SD_A, SD_E, Scale, Range_X, Range_Y, SpatialSimModel, n_s, Ngrid_sim, M, 
  RecFn, alpha_g, ro, MRPSB, F_equil, S_bioecon, Accel, SD_F, n_samp_per_year, AreaSwept, q_I, mu_R0_total){

  # Derived values
  DomainArea = diff(Range_X) * diff(Range_Y)
  w_a = rep(NA,length=20)
    w_a[1] = alpha_g
    for(a in 2:length(w_a)) w_a[a] = alpha_g + ro*w_a[a-1]
    w_k = w_a[k]
    MRPSB = (1 - exp(-M) - ro*exp(-M) + ro*exp(-2*M) ) / ( w_k - (w_k-alpha_g)*exp(-M) )  # Maximum recruits per spawning biomass
  mu_R0_per_area = mu_R0_total / DomainArea
    mu_S0_per_area = mu_R0_per_area / MRPSB
  if(RecFn=="Ricker"){
    log_mu_alpha = log(mu_R0_total)   # median of maximum recruits per spawning biomass per unit area; min possible = MRPSB
    beta = log( mu_R0/exp(log_mu_alpha)/mu_S0_per_area ) / (-mu_S0_per_area)
  }
  if(RecFn=="BH"){
    log_mu_alpha = log(mu_R0_total)   # median of maximum recruits per spawning biomass per unit area; min possible = MRPSB
    beta = (exp(log_mu_alpha)*mu_S0_per_area/mu_R0_per_area - 1) / mu_S0_per_area
    #tmp = seq(0,mu_S0_per_area,length=1000)
    #tmp2 = exp(log_mu_alpha)*tmp / (1 + beta*tmp)
    #plot( x=tmp, y=tmp2) 
  }
  if(RecFn=="Constant"){
    log_mu_alpha = log(mu_R0_total / DomainArea) # Expected recruitment per unit area
  }

  # Simulate random variability
  Eps_F = rnorm(n_t, mean=0, sd=SD_F)
  Eps_C = rnorm(n_t, mean=0, sd=CV_c)
  Year_Range = c(1,n_t)
  F_t = rep(NA, n_t)

  # Define fields
  Nu=1
  if(SpatialSimModel=="Matern"){
    # matern: C(r=h/Scale) = Var * 2^{1-Nu} * gamma(Nu)^{-1} * (sqrt{2*Nu}*r)^Nu * besselK(sqrt{2*Nu} * r)
    model_A <- RMmatern(nu=Nu, var=SD_A^2, scale=Scale)
    model_E <- RMmatern(nu=Nu, var=SD_E^2, scale=Scale)
  }
  # Calculate range at which correlation is 13% of maximum (comparable to range defined in Lindgren and Rue 2013, immediately before Eq. 4)
  h = seq(-3*diff(Range_X),3*diff(Range_X),length=10000)
  r = abs(h)/Scale
  if(SpatialSimModel=="Matern") C = SD_A^2 * 2^(1-Nu) * gamma(Nu)^(-1) * (sqrt(2*Nu)*r)^Nu * besselK(sqrt(2*Nu) * r, nu=Nu)
  if(SpatialSimModel=="Gaussian") C = SD_A^2 * exp(-r^2)
  if(FALSE){
    par(mfrow=c(1,2))
    plot(model_A, ylim=c(0,1), xlim=c(-3,3)*diff(Range_X))
    plot(x=h, y=C, ylim=c(0,1), type="l")
  }
  RangeTrue = abs(h[which.min( (C-0.10*SD_A)^2 )])

  
  # locations of sampling stations
  if( n_samp_per_year < n_s*2 ) error("Increase n_samp_per_year for K-means algorithm")
  if(TRUE){
    x_stations = runif(n_samp_per_year, min=Range_X[1], max=Range_X[2])
    y_stations = runif(n_samp_per_year, min=Range_Y[1], max=Range_Y[2])
    loc_samples = cbind("long"=x_stations, "lat"=y_stations)
  }
  if(FALSE){
    x_stations = seq(Range_X[1],Range_X[2],length=sqrt(n_samp_per_year))
    y_stations = seq(Range_Y[1],Range_Y[2],length=sqrt(n_samp_per_year))
    loc_samples = expand.grid("long"=x_stations, "lat"=y_stations)
  }

  # K-means
  K = kmeans( x=loc_samples[,c('long','lat')], centers=n_s, iter.max=1e3, nstart=100) # K$tot.withinss
    K$centers = K$centers + rnorm( n_s, mean=0, sd=rep(0.001*c(diff(Range_X),diff(Range_Y)),each=n_s) )  # Necessary to ensure Kmeans aren't identical to other samples

  # Generate list of points to track
  loc_samples = cbind( loc_samples, "Station_j"=K$cluster )
  loc_stations = cbind( K$centers, "Station_j"=1:n_s )
  loc_grid = expand.grid( "long"=seq(Range_X[1],Range_X[2],length=ceiling(sqrt(Ngrid_sim))),"lat"=seq(Range_Y[1],Range_Y[2],length=ceiling(sqrt(Ngrid_sim))) )
    NN = nn2( query=loc_grid, data=loc_stations[,1:2], k=1 )
    loc_grid = cbind( loc_grid, "Station_j"=NN$nn.idx )
  loc_all = rbind( loc_stations, loc_samples, loc_grid )

  # Get area for each location
  Voronoi_samples = calcVoronoi( cbind(X=loc_samples[,'long'], Y=loc_samples[,'lat']), xlim=range(Range_X,loc_samples[,'long']), ylim=range(Range_Y,loc_samples[,'lat']))
  Area_samples = calcArea( Voronoi_samples )[,2]
  Area_all = c( rep(0,n_s), Area_samples, rep(0,Ngrid_sim) )

  # Simulate null field
  Omega_s = RFsimulate(model=model_A, x=loc_all[,'long'], y=loc_all[,'lat'])@data[,1]
    Alpha_s = exp(Omega_s + log_mu_alpha)

  # Annual variation
  Epsilon_tmp_s = Epsilon_s = array(NA, dim=c(n_s+n_samp_per_year+Ngrid_sim, n_t))
  for(t in 1:n_t){
    Epsilon_s[,t] = Epsilon_tmp_s[,t] = RFsimulate(model_E, x=loc_all[,'long'], y=loc_all[,'lat'])@data[,1]
  }

  # Equilibrium conditions - recruitment
  if(RecFn=="Ricker"){
    mu_S_equil_per_area = log((1 - exp(-M-F_equil) - ro*exp(-M-F_equil) + ro*exp(-2*M-2*F_equil) ) / (exp(log_mu_alpha) * ( w_k - (w_k-alpha_g)*exp(-M-F_equil) ))) / (-beta)
    mu_R_equil_per_area = exp(log_mu_alpha)*mu_S_equil*exp(-beta*mu_S_equil)
  }
  if(RecFn=="BH"){
    mu_S_equil_per_area = ((exp(log_mu_alpha) * ( w_k - (w_k-alpha_g)*exp(-M-F_equil) ) / (1 - exp(-M-F_equil) - ro*exp(-M-F_equil) + ro*exp(-2*M-2*F_equil) )) - 1) / beta
    mu_R_equil_per_area = exp(log_mu_alpha)*mu_S_equil / (1+beta*mu_S_equil)
  }
  if(RecFn=="Constant"){
    mu_S_equil_per_area = (exp(log_mu_alpha) * ( w_k - (w_k-alpha_g)*exp(-M-F_equil) ) / (1 - exp(-M-F_equil) - ro*exp(-M-F_equil) + ro*exp(-2*M-2*F_equil) ))
    mu_R_equil_per_area = exp(log_mu_alpha)
    mu_S0_per_area = (exp(log_mu_alpha) * ( w_k - (w_k-alpha_g)*exp(-M-(0)) ) / (1 - exp(-M-(0)) - ro*exp(-M-(0)) + ro*exp(-2*M-2*(0)) ))
  }
  mu_N_equil_per_area = mu_R_equil_per_area / (1 - exp(-M-F_equil))
  S_equil = mu_S_equil_per_area * exp( Omega_s )
  R_equil = mu_R_equil_per_area * exp( Omega_s )
  N_equil = mu_N_equil_per_area * exp( Omega_s )
  W_equil = S_equil / N_equil
  S0 = mu_S0_per_area * exp( Omega_s )
  sum_S0 = sum(S0 * Area_all)
  if( abs(log(sum_S0 / (sum( R_equil * Area_all )/MRPSB)))>0.5  ) stop("Problem with area units")
  print( paste("Expected unfished catch rate = ", sum( mu_S_equil_per_area*AreaSwept ) ))

  #### Simulate B
  C_st = S_st = N_st = R_st = W_st = matrix(NA, nrow=n_s+n_samp_per_year+Ngrid_sim, ncol=n_t) # These are "per unit area"
  # Project forward
  for(t in 1:n_t){
    # Recruitment
    R_st[,t] = R_equil * exp(Epsilon_s[,t])
    # Numbers
    if(t==1) N_st[,t] = N_equil * exp(-M-F_equil) + R_st[,t]
    if(t>=2) N_st[,t] = N_st[,t-1] * exp(-M-F_t[t-1]) + R_st[,t]
    # Biomass
    #if(t==1) S_st[,t] = (1+ro)*exp(-M-F_equil)*S_equil - ro*exp(-2*M-F_equil-F_equil)*S_equil + w_k*R_st[,t] - (w_k-alpha_g)*exp(-F_equil-M)*R_equil
    #if(t==2) S_st[,t] = (1+ro)*exp(-M-F_t[t-1])*S_st[,t-1] - ro*exp(-2*M-F_t[t-1]-F_equil)*S_equil + w_k*R_st[,t] - (w_k-alpha_g)*exp(-F_t[t-1]-M)*R_st[,t-1]
    #if(t>=3) S_st[,t] = (1+ro)*exp(-M-F_t[t-1])*S_st[,t-1] - ro*exp(-2*M-F_t[t-1]-F_t[t-2])*S_st[,t-2] + w_k*R_st[,t] - (w_k-alpha_g)*exp(-F_t[t-1]-M)*R_st[,t-1]
    if(t==1) S_st[,t] = exp(-M-F_equil)* (alpha_g*N_equil + ro*S_equil) + w_k*R_st[,t]
    if(t>=2) S_st[,t] = exp(-M-F_t[t-1])* (alpha_g*N_st[,t-1] + ro*S_st[,t-1]) + w_k*R_st[,t]
    # Weight
    W_st[,t] = S_st[,t] / N_st[,t]
    # Fishing mortality
    if(t==1) F_t[t] = F_equil
    if(t>=2) F_t[t] = F_t[t-1] * (sum(Area_all*S_st[,t])/sum(Area_all*S_equil)/S_bioecon)^Accel * exp( Eps_F[t] - SD_F^2/2 )
    # Catches
    C_st[,t] = F_t[t]/(M+F_t[t]) * (1-exp(-M-F_t[t])) * S_st[,t]
  }
  ## Simulate data
  # Survey catches
  Y = I_st = rpois(n=n_samp_per_year*n_t, lambda=AreaSwept*q_I*N_st[n_s+1:n_samp_per_year,] )
  Y_stations = matrix(Y, ncol=n_t, byrow=FALSE)
  NAind_stations = array( as.integer(ifelse(is.na(Y_stations),1,0)), dim=dim(Y_stations))
  # Average weight
  Wobs_st = array(rnorm(n=n_samp_per_year*n_t, mean=W_st[n_s+1:n_samp_per_year,], sd=W_st[n_s+1:n_samp_per_year,]*CV_w), dim=c(n_samp_per_year,n_t))
  # Total catch
  C_t = colSums(C_st * outer(Area_all,rep(1,ncol(C_st))))
  Cobs_t = (C_t + C_t*Eps_C)
  # Total abundance
  S_t = colSums(S_st*outer(Area_all,rep(1,ncol(C_st))))
  R_t = colSums(R_st*outer(Area_all,rep(1,ncol(C_st))))
  N_t = colSums(N_st*outer(Area_all,rep(1,ncol(C_st))))

  # Convert data to data frame
  DF = data.frame( "I_j"=as.vector(I_st), "W_j"=as.vector(Wobs_st), "Station_j"=rep(loc_samples[,'Station_j'],n_t), "Year_j"=as.vector(col(Wobs_st)), "AreaSwept_j"=AreaSwept, "long_j"=rep(loc_samples[,'long'],n_t), "lat_j"=rep(loc_samples[,'lat'],n_t) )

  # List of generated data
  Return = list("Eps_F"=Eps_F, "Eps_C"=Eps_C, "Year_Range"=Year_Range, "F_t"=F_t, "RangeTrue"=RangeTrue, 
  "loc_samples"=loc_samples, "K"=K, "loc_samples"=loc_samples, "loc_stations"=loc_stations, "loc_grid"=loc_grid, 
  "NN"=NN, "loc_grid"=loc_grid, "loc_all"=loc_all, "Area_samples"=Area_samples, "Area_all"=Area_all, 
  "Voronoi_samples"=Voronoi_samples, "Omega_s"=Omega_s, "Alpha_s"=Alpha_s, "Epsilon_s"=Epsilon_s, 
  "mu_S_equil_per_area"=mu_S_equil_per_area, "mu_R_equil_per_area"=mu_R_equil_per_area, "mu_S0_per_area"=mu_S0_per_area, 
  "mu_N_equil_per_area"=mu_N_equil_per_area, "S_equil"=S_equil, "R_equil"=R_equil, "N_equil"=N_equil, 
  "W_equil"=W_equil, "S0"=S0, "sum_S0"=sum_S0, "C_st"=C_st, "S_st"=S_st, "N_st"=N_st, "R_st"=R_st, "W_st"=W_st, 
  "I_st"=I_st, "Y"=Y, "Y_stations"=Y_stations, "NAind_stations"=NAind_stations, "Wobs_st"=Wobs_st, "C_t"=C_t, 
  "Cobs_t"=Cobs_t, "S_t"=S_t, "R_t"=R_t, "N_t"=N_t, "w_a"=w_a, "w_k"=w_k, "MRPSB"=MRPSB, 
  "mu_R0_per_area"=mu_R0_per_area, "mu_S0_per_area"=mu_S0_per_area, "log_mu_alpha"=log_mu_alpha, "beta"=beta,
  "DomainArea"=DomainArea, "DF"=DF)
  return( Return )

}
