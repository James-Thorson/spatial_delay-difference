MakeTmbList_Fn <-
function( Version, Model, Fix_Q, ErrorModel_CatchRates, ErrorModel, Smooth_F, n_j, n_i, n_s, n_t, DF_input, 
  C_t, IndexMat, mesh_stations, spde_stations, Area_i, alpha_g, ro, w_k, M, k, CV_c, CV_w, q_I ){

  # Starting values
  F_t = rep(0.1, n_t)

  # Version-specific stuff
  if(Version=="delay_difference_v4d"){
    Error_Model = ErrorModel_CatchRates
    Data = list(Error_Model=Error_Model, Smooth_F=Smooth_F, n_j=n_j, n_i=n_i, n_s=n_s, n_t=n_t, I_j=DF_input[,'I_j'], W_j=DF_input[,'W_j'], AreaSwept_j=DF_input[,'AreaSwept_j'], Station_j=DF_input[,'Station_j']-1, Year_j=DF_input[,'Year_j']-1, C_t=C_t, meshidxloc=mesh_stations$idx$loc-1, G0=spde_stations$param.inla$M0, G1=spde_stations$param.inla$M1, G2=spde_stations$param.inla$M2, Area_i=Area_i, alpha_g=alpha_g, ro=ro, w_k=w_k, M=M, k=k, CV_c=CV_c, CV_w=CV_w)
    Parameters = list(log_F_sd=log(1), log_F_equil=log(F_t[1]), log_F_t_input=log(F_t), log_q_I=log(q_I), beta=c(0.0), log_tau_E=0.0, log_tau_O=0.0, log_kappa=0.0,	ln_VarInfl=c(0.0,0.0), Epsilon_input=matrix(0,spde_stations$n.spde,n_t), Omega_input=rep(0,spde_stations$n.spde))
  }
  if(Version=="delay_difference_v6c"){
    Data = list(ModelType=switch(Model,"Spatial"=1,"Nonspatial"=2,"Index"=2,"MeanLength_terminalcatch"=2), ErrorModel_CatchRates=ErrorModel_CatchRates, ErrorModel_MeanWeight=ErrorModel_MeanWeight, Smooth_F=Smooth_F, n_j=n_j, n_i=n_i, n_s=n_s, n_t=n_t, I_j=DF_input[,'I_j'], W_j=DF_input[,'W_j'], AreaSwept_j=DF_input[,'AreaSwept_j'], Station_j=DF_input[,'Station_j']-1, Year_j=DF_input[,'Year_j']-1, C_t=C_t, meshidxloc=mesh_stations$idx$loc-1, G0=spde_stations$param.inla$M0, G1=spde_stations$param.inla$M1, G2=spde_stations$param.inla$M2, Area_i=Area_i, alpha_g=alpha_g, ro=ro, w_k=w_k, M=M, k=k, CV_c=CV_c, CV_w=CV_w)
    Parameters = list(log_F_sd=log(1), log_F_equil=log(F_equil), log_F_t_input=log(F_t), log_q_I=log(q_I), beta=c(0.0), log_tau_E=0.0, log_tau_O=0.0, log_kappa=0.0,	ln_VarInfl=c(0.0,0.0), log_extraCV_w=log(0.05), log_tau_N=log(1), Epsilon_input=matrix(0,spde_stations$n.spde,n_t), Omega_input=rep(0,spde_stations$n.spde), Nu_input=rep(0,n_t))
  }
  if(Version %in% c("delay_difference_v8e","delay_difference_v8d")){
    Data = list(ModelType=switch(Model,"Spatial"=1,"Nonspatial"=2,"Index"=2,"MeanLength_terminalcatch"=2), ErrorModel_CatchRates=ErrorModel_CatchRates, ErrorModel_MeanWeight=ErrorModel_MeanWeight, Smooth_F=Smooth_F, n_j=n_j, n_i=n_i, n_s=n_s, n_t=n_t, I_j=DF_input[,'I_j'], W_j=DF_input[,'W_j'], AreaSwept_j=DF_input[,'AreaSwept_j'], Station_j=DF_input[,'Station_j']-1, Year_j=DF_input[,'Year_j']-1, C_t=C_t, IndexMat=IndexMat, meshidxloc=mesh_stations$idx$loc-1, G0=spde_stations$param.inla$M0, G1=spde_stations$param.inla$M1, G2=spde_stations$param.inla$M2, Area_i=Area_i, alpha_g=alpha_g, ro=ro, w_k=w_k, M=M, k=k, CV_c=CV_c, CV_w=CV_w)
    Parameters = list(log_F_sd=log(1), log_F_t_input=log(c(F_equil,F_t)), log_q_I=log(q_I), beta=log_mu_alpha, log_tau_E=0.0, log_tau_O=0.0, log_kappa=0.0,	ln_VarInfl=c(0.0,0.0), log_extraCV_w=log(0.05), log_tau_N=log(1), log_extraCV_Index=rep(log(0.1),2), Epsilon_input=matrix(0,spde_stations$n.spde,n_t), Omega_input=rep(0,spde_stations$n.spde), Nu_input=rep(0,n_t))
  }

  # Define random effects
  if( Model=="Spatial" ){
    if( Smooth_F==0 ) Random = c("Epsilon_input", "Omega_input")
    if( Smooth_F!=0 ) Random = c("Epsilon_input", "Omega_input", "log_F_equil")
  }
  if( Model=="Nonspatial" ){
    if( Smooth_F==0 ) Random = c("Nu_input")
    if( Smooth_F!=0 ) Random = c("Nu_input", "log_F_t_input")
  }
  if( Model=="Index" ){
    if( Smooth_F==0 ) Random = c("Nu_input")
    if( Smooth_F!=0 ) Random = c("Nu_input", "log_F_t_input")
  }
  if( Model=="MeanLength_terminalcatch" ){
    if( Smooth_F==0 ) Random = c("Nu_input")
    if( Smooth_F!=0 ) Random = c("Nu_input", "log_F_t_input")
  }

  # Define fixed parameters
  Map = list()
  if( Smooth_F==0 ) Map[["log_F_sd"]] = factor(NA)
  if( ErrorModel_CatchRates==0 ) Map[["ln_VarInfl"]] = factor(c(NA,NA))
  if( ErrorModel_MeanWeight==0 ) Map[["log_extraCV_w"]] = factor(NA)
  if( Model=="Spatial" ){
    Map[["log_tau_N"]] = factor(NA)
    Map[["Nu_input"]] = factor( rep(NA,n_t) )
    Map[["log_extraCV_Index"]] = factor( rep(NA,2) )
    Map[["log_F_t_input"]] = factor( c(1,1:n_t) )
  }
  if( Model=="Nonspatial" ){
    Map[["log_tau_E"]] = factor(NA)
    Map[["log_tau_O"]] = factor(NA)
    Map[["Epsilon_input"]] = factor( rep(NA,length(Parameters$Epsilon_input)) )
    Map[["Omega_input"]] = factor( rep(NA,length(Parameters$Omega_input)) )
    Map[["log_kappa"]] = factor(NA)
    Map[["log_extraCV_Index"]] = factor( rep(NA,2) )
    Map[["log_F_t_input"]] = factor( c(1,1:n_t) )
  }
  if( Model=="Index" ){
    Map[["log_tau_E"]] = factor(NA)
    Map[["log_tau_O"]] = factor(NA)
    Map[["Epsilon_input"]] = factor( rep(NA,length(Parameters$Epsilon_input)) )
    Map[["Omega_input"]] = factor( rep(NA,length(Parameters$Omega_input)) )
    Map[["log_kappa"]] = factor(NA)
    Map[["ln_VarInfl"]] = factor(c(NA,NA))
    Map[["log_extraCV_w"]] = factor(NA)
    Map[["log_F_t_input"]] = factor( c(1,1:n_t) )
  }
  if( Model=="MeanLength_terminalcatch" ){
    Map[["log_tau_E"]] = factor(NA)
    Map[["log_tau_O"]] = factor(NA)
    Map[["Epsilon_input"]] = factor( rep(NA,length(Parameters[["Epsilon_input"]])) )
    Map[["Omega_input"]] = factor( rep(NA,length(Parameters[["Omega_input"]])) )
    Map[["log_kappa"]] = factor(NA)
    Map[["ln_VarInfl"]] = factor(c(NA,NA))
    Map[["log_extraCV_w"]] = factor(NA)
    Map[["log_F_t_input"]] = factor( c(1,1:n_t) )
    # Turn off variance inflation
    Map[["log_extraCV_Index"]] = factor( c(NA,NA) )
    Parameters[["log_extraCV_Index"]] = rep(-5,2)
    # Turn off rec devs
    Map[["log_tau_N"]] = factor(NA)
    Map[["Nu_input"]] = factor( rep(NA,length(Parameters[["Nu_input"]])) )
    # Fix variation in F
    Map[["log_F_sd"]] = factor(NA)
    Parameters[["log_F_sd"]] = 0.1 
  }
  if( Fix_Q==TRUE ) Map[["log_q_I"]] = factor(NA)
  
  # Exclude data depending upon which model setting to use
  if(Model == "Index"){
    Data[["I_j"]][] = -999
    Data[["W_j"]][] = -999
  }
  if(Model %in% c("Nonspatial","Spatial")){
    Data[["IndexMat"]][] = -999
  }
  if(Model == "MeanLength_terminalcatch"){
    Data[["I_j"]][] = -999
    Data[["W_j"]][] = -999
    Data[["C_t"]][-Data$n_t] = -999
    Data[["IndexMat"]][,c('Mean_I','logSD_I')] = -999
  }
  
  # Return objects
  Return = list("Parameters"=Parameters, "Data"=Data, "Map"=Map, "Random"=Random)
  return( Return )
}
