

#########################
# Simulation example
#########################

# Install from GitHub
library(devtools)
install_github("James-Thorson/spatial_delay-difference")

RunFile = ### Set this to something

# Location of TMB files
TmbFile = paste0(system.file("executables", package="SpatialDelayDiff"),"/")

# Libraries
library(INLA)
library(TMB)
  newtonOption(smartsearch=TRUE)
library(SpatialDelayDiff)

# Which model to use?
  Model = c("Nonspatial", "Spatial", "Index", "MeanLength_terminalcatch")[2]

# Compile model in TMB
  setwd( TmbFile )
  Version = "delay_difference_v8e"
  compile( paste0(Version,".cpp") )

#### Settings
# Simulation test
  RepSet = 1:200
  RepSet = RepSet + max(RepSet)*(ThreadNum-1) 
# Estimation
  ErrorModel_CatchRates = 1   # 0: Poisson; 1: Negative binomial for counts
  ErrorModel_MeanWeight = 1   # 0: Fixed-CV; 1: Est-CV
  Smooth_F = 0      # 0: No; 1: Yes
  Fix_Q = TRUE   # 
  SpatialSimModel = "Matern"
  MeshType = c("Minimal","Recommended")[1]
  n_s = 25
# Domain
  n_t = 30
  Range_X = c(0,1000)
  Range_Y = c(0,1000)
# Biological - growth
  k = 3         # Age at when individual is fully recruited
  ro = 0.8      # Brody growth coefficient
  alpha_g = 0.2   # Weight at age 1
# Biological - survival
  M = 0.3  # Natural mortality 
# Biological recruitment
  RecFn = c("Ricker", "BH", "Constant")[3]
  mu_R0_total = 1e9 # Median total R0
    mu_N0 = mu_R0_total / (1-exp(-M))
# Variability
  SD_A = 0.5  # 0.5  # Spatial variation in productivity
  SD_E = 0.5  # 0.5  # Spatiotemporal variation in recruitment
  Scale = 0.25 * 1000
# Fishing mortality 
  F_equil = 0.05  # Initial equilibrium fishing mortality
  S_bioecon = 0.4
  Accel = 0.2
  SD_F = 0.2 # 0.2
# Data
  n_samp_per_year = 100 
  AreaSwept = 0.0025 # 10 / mu_N0 * DomainArea  # in km^2
  q_I = 1
  CV_w = 0.20
  CV_c = 0.05
# Visualization
  Ngrid_sim = 1e4
  Ngrid_proj = 1e4  
  
  # Save settings
  SettingsList = list( "RepSet"=RepSet, "ErrorModel_CatchRates"=ErrorModel_CatchRates, "ErrorModel_MeanWeight"=ErrorModel_MeanWeight, "Smooth_F"=Smooth_F, "Fix_Q"=Fix_Q, "SpatialSimModel"=SpatialSimModel, "MeshType"=MeshType, "Version"=Version, "n_s"=n_s, "n_t"=n_t, "Range_X"=Range_X, "Range_Y"=Range_Y, "k"=k, "ro"=ro, "alpha_g"=alpha_g, "M"=M, "RecFn"=RecFn, "mu_R0_total"=mu_R0_total, "SD_A"=SD_A, "SD_E"=SD_E, "Scale"=Scale, "F_equil"=F_equil, "S_bioecon"=S_bioecon, "Accel"=Accel, "SD_F"=SD_F, "n_samp_per_year"=n_samp_per_year, "AreaSwept"=AreaSwept, "q_I"=q_I, "CV_w"=CV_w, "CV_c"=CV_c)
    capture.output(SettingsList, file=paste(RunFile,"SettingsList.txt",sep=""))
    save(SettingsList, file=paste(RunFile,"SettingsList.RData",sep=""))
    file.copy( from=paste(TmbFile,Version,".cpp",sep=""), to=paste(RunFile,Version,".cpp",sep=""), overwrite=TRUE)
  
  # Simulate data
  DataList = SimData_Fn(n_t=n_t, CV_C=CV_C, SD_A=SD_A, SD_E=SD_E, Scale=Scale, Range_X=Range_X, 
    Range_Y=Range_Y, SpatialSimModel=SpatialSimModel, n_s=n_s, Ngrid_sim=Ngrid_sim,
    M=M, RecFn=RecFn, alpha_g=alpha_g, ro=ro, MRPSB=MRPSB, F_equil=F_equil, S_bioecon=S_bioecon, 
    Accel=Accel, SD_F=SD_F, n_samp_per_year=n_samp_per_year, AreaSwept=AreaSwept, q_I=q_I, mu_R0_total=mu_R0_total)
  # Attach data
  attach(DataList)

  # Make INLA inputs for TMB
  if(MeshType=="Recommended"){
    mesh_stations = inla.mesh.create( loc_stations[,c('lat','long')], plot.delay=NULL, extend=list(n=8,offset=-0.15), refine=list(min.angle=26,max.edge.data=0.08,max.edge.extra=0.2) )  # loc_samp
  }
  if(MeshType=="Minimal"){
    mesh_stations = inla.mesh.create( loc_stations[,c('lat','long')], plot.delay=NULL, extend=list(n=8,offset=-0.15), refine=F )  # loc_samp
  }
  spde_stations = inla.spde2.matern(mesh_stations, alpha=2)  # Given 2D field, alpha=2 <=> Matern Nu=1 (Lindgren and Rue 2013, between Eq. 1&2)
  n_i = nrow(mesh_stations$loc)
  n_j = nrow(DF)

  # Get area for each location
  Voronoi_s = calcVoronoi( cbind(X=loc_stations[,1], Y=loc_stations[,2]), xlim=range(Range_X,loc_stations[,1]), ylim=range(Range_Y,loc_stations[,2]))
  Area_s = calcArea( Voronoi_s )[,2]
  Area_i = c( Area_s, rep(0,n_i-n_s) )

  # Calculate index of abundance
  PredDF = data.frame( "Year_j"=sort(unique(DF[,'Year_j'])) )
    Glm_I = zeroinfl(I_j ~ 0 + factor(Year_j) | 1, data=DF, dist=c("poisson","negbin","geometric")[2], method="Nelder-Mead", control=zeroinfl.control("EM"=TRUE))
    Glm_W = lm(W_j ~ 0 + factor(Year_j), data=DF, subset=DF[,'I_j']>0)
    Index_hat = array(-999, dim=c(n_t,4), dimnames=list(Year_Range[1]:Year_Range[2],c("Mean_I","logSD_I","Mean_W","SD_W")))
    Index_hat[unique(DF[,'Year_j']),c("Mean_I","logSD_I")] = cbind(predict(Glm_I, newdata=PredDF, type="response")/(AreaSwept*q_I)*DomainArea, summary(Glm_I)$coef$count[1:length(unique(DF[,'Year_j'])),'Std. Error'])
    Index_hat[unique(DF[,'Year_j']),c("Mean_W","SD_W")] = cbind(predict(Glm_W, newdata=PredDF, type="response"), summary(Glm_W)$coef[1:length(unique(DF[,'Year_j'])),'Std. Error'])
  
  # Make TMB inputs
  TmbList = MakeTmbList_Fn( Version=Version, Model=Model, IndexMat=Index_hat, Fix_Q=Fix_Q, ErrorModel_CatchRates=ErrorModel_CatchRates, ErrorModel=ErrorModel, Smooth_F=Smooth_F, n_j=n_j, n_i=n_i, n_s=n_s, n_t=n_t, DF_input=DF, C_t=C_t, mesh_stations=mesh_stations, spde_stations=spde_stations, Area_i=Area_i, alpha_g=alpha_g, ro=ro, w_k=w_k, M=M, k=k, CV_c=CV_c, CV_w=CV_w, q_I=q_I )
  # Look at data that are passed for this model
  # TmbList$Data[-match(c("G0","G1","G2"),names(TmbList$Data))]
  
  # Save stuff
  MapList = list( "Map"=TmbList[["Map"]], "Random"=TmbList[["Random"]])
  capture.output(MapList, file=paste(RunFile,"MapList.txt",sep=""))

  # Build object                                                              #    
  dyn.load( paste(TmbFile,dynlib(Version),sep="") )
  obj <- MakeADFun(data=TmbList[["Data"]], parameters=TmbList[["Parameters"]], random=TmbList[["Random"]], map=TmbList[["Map"]], hessian=FALSE)
  
  # First runs
  obj$fn( obj$par )
  obj$gr( obj$par )
  
  # Look at which parameters are being estimated
  obj$env$last.par.best
  
  # Settings
  obj$control <- list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=100)
  obj$hessian <- FALSE

  # Run optimizer
  opt = nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, lower=-20, upper=20, control=list(trace=1, eval.max=1e4, iter.max=1e4))
  opt[["final_gradient"]] = obj$gr( opt$par )
  Report = obj$report()
  Sdreport = try( sdreport(obj) )
  dyn.unload( paste(TmbFile,dynlib(Version),sep="") )
  
  # Save results
  capture.output( opt, file=paste(RunFile,"opt.txt",sep=""))
  capture.output( Report, file=paste(RunFile,"Report.txt",sep=""))
  Save = list( "opt"=opt, "Report"=Report, "obj"=obj, "Sdreport"=Sdreport)
  save( Save, file=paste(RunFile,"Save.RData",sep=""))

  ############################
  # Plots to visualize fit
  ############################

  # Save Varanoi stuff
  png( file=paste(RunFile,"Voronoi.png",sep=""), width=6, height=3, res=200, units="in")
    par( mfrow=c(1,2), mar=c(3,3,2,0), mgp=c(2,0.5,0) )
    # Visualize samples
    plotMap( Voronoi_samples )
    points( y=loc_samples[,'lat'], x=loc_samples[,'long'], col=c("red","black")[rep.int( c(1,2), times=c(n_s,n_samp_per_year-n_s))], pch=loc_samples[,'Station_j'] )
    # Visualize stations
    plotMap( Voronoi_s )
    points( y=loc_samples[,'lat'], x=loc_samples[,'long'], col=c("red","black")[rep.int( c(1,2), times=c(n_s,n_samp_per_year-n_s))] )
  dev.off()

  # Plot time series
  png( file=paste(RunFile,"True_Timeseries.png",sep=""), width=3*3, height=2*3, res=200, units="in")
    par(mfrow=c(2,3), mar=c(3,3,2,0))
    plot( C_t, type="l", ylim=c(0,max(C_t)), main="Catch" )
    plot( F_t, type="l", ylim=c(0,max(F_t)), main="Fishing mortality" )
    plot( S_t, type="l", ylim=c(0,max(S_t)), main="Spawning biomass" )
    plot( S_t/sum_S0, type="l", ylim=c(0,1.2), main="Depletion" )
    plot( N_t, type="l", ylim=c(0,max(N_t)), main="Total abundance" )
    for(t in 1:n_t){
      points( x=t, y=Index_hat[t,'Mean_I'], col="red")
      lines(x=c(t,t), y=Index_hat[t,'Mean_I']*exp(c(-1,1)*Index_hat[t,'logSD_I']), col="red" )
    }
    plot( S_t/N_t, type="l", ylim=c(0,max(S_t/N_t)), main="Average weight" )
    for(t in 1:n_t){
      points( x=t, y=Index_hat[t,'Mean_W'], col="red")
      lines(x=c(t,t), y=Index_hat[t,'Mean_W']+c(-1,1)*Index_hat[t,'SD_W'], col="red" )
    }
  dev.off()
  
  # Check correlation of random fields
  Cor_Omega = cor( Report$Omega[1:n_s], Omega_s[1:n_s] )
  Cor_Epsilon = rep(NA, n_t)
  for(t in 1:n_t) Cor_Epsilon[t] = cor( Report$Epsilon[1:n_s,t], Epsilon_s[1:n_s,t] )
  
  # Colors
  Col = colorRampPalette(colors=c("blue","grey","red"))
  f = function(Num){ if( var(as.vector(Num),na.rm=TRUE)==0 ){ Return = array(rep(0.5,prod(dim(Num))),dim(Num))}else{Return=(plogis(Num)-min(plogis(Num),na.rm=TRUE))/diff(range(plogis(Num),na.rm=TRUE))}; return(Return) } 
  RangeFn = function(Vec) range(ifelse(abs(Vec)==Inf,NA,Vec),na.rm=TRUE) 
        
  # Time series
  FUN = function(InputMat) c( InputMat[,1]-1*InputMat[,2], rev(InputMat[,1]+1*InputMat[,2]))
  png( file=paste(RunFile,"True_vs_Est.png",sep=""), width=9, height=6, res=200, units="in")
    par(mfrow=c(2,3), mar=c(3,3,2,0))
    # Abundance
    Mat = cbind( "Year"=unique(DF[,'Year_j']), "True"=N_t, "Est"=Report$N_t_hat)
    matplot( y=Mat[,c("True",'Est')], x=Mat[,c("Year")], type="l", col=c("black","red"), lty="solid", ylim=c(0,max(Mat[,c('True','Est')])), main="Abundance")
    if( !("condition" %in% names(attributes(Sdreport))) ) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="N_t_hat"),]), x=c(Mat[,c("Year")],rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
    # Biomass
    Mat = cbind( "Year"=unique(DF[,'Year_j']), "True"=S_t, "Est"=Report$S_t_hat)
    matplot( y=Mat[,c("True",'Est')], x=Mat[,c("Year")], type="l", col=c("black","red"), lty="solid", ylim=c(0,max(Mat[,c('True','Est')])), main="Spawning biomass")
    if( !("condition" %in% names(attributes(Sdreport))) ) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_t_hat"),]), x=c(Mat[,c("Year")],rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
    # Average weight
    Mat = cbind( "Year"=unique(DF[,'Year_j']), "True"=S_t/N_t, "Est"=Report$S_t_hat/Report$N_t_hat)
    matplot( y=Mat[,c("True",'Est')], x=Mat[,c("Year")], type="l", col=c("black","red"), lty="solid", ylim=c(0,max(Mat[,c('True','Est')])), main="Average weight")
    if( !("condition" %in% names(attributes(Sdreport))) ) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_t_hat / N_t_hat"),]), x=c(Mat[,c("Year")],rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
    # Recruitment
    Mat = cbind( "Year"=unique(DF[,'Year_j']), "True"=R_t, "Est"=Report$R_t_hat)
    matplot( y=Mat[,c("True",'Est')], x=Mat[,c("Year")], type="l", col=c("black","red"), lty="solid", ylim=c(0,max(Mat[,c('True','Est')])), main="Recruitment")
    if( !("condition" %in% names(attributes(Sdreport))) ) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="R_t_hat"),]), x=c(Mat[,c("Year")],rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
    # Fishing mortality
    Mat = cbind( "Year"=unique(DF[,'Year_j']), "True"=F_t, "Est"=Report$F_t)
    matplot( y=Mat[,c("True",'Est')], x=Mat[,c("Year")], type="l", col=c("black","red"), lty="solid", ylim=c(0,max(Mat[,c('True','Est')])), main="Fishing mortality")
    if( !("condition" %in% names(attributes(Sdreport))) ) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="F_t"),]), x=c(Mat[,c("Year")],rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
    # Depletion
    Mat = cbind( "Year"=unique(DF[,'Year_j']), "True"=S_t/sum_S0, "Est"=Report$S_t/Report$sum_S0)
    matplot( y=Mat[,c("True",'Est')], x=Mat[,c("Year")], type="l", col=c("black","red"), lty="solid", ylim=c(0,max(Mat[,c('True','Est')])), main="Depletion")
    if( !("condition" %in% names(attributes(Sdreport))) ) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="S_t_hat / sum_S0"),]), x=c(Mat[,c("Year")],rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)  
  dev.off()

  # fields
  Years2Show = 1:3
  png( file=paste(RunFile,"True_vs_Est_Recruitment.png",sep=""), width=2*2, height=(1+length(Years2Show))*2, res=200, units="in")
    par(mfrow=c(1+length(Years2Show),2), mar=c(0.2,0.2,0,0), oma=c(4,4,2,2), mgp=c(2,0.5,0), tck=-0.02)
    Zlim = range( log(c(R_st[-c(1:(n_s+n_samp_per_year)),Years2Show],Report$R_it[1:n_s,Years2Show])) ) 
    X = seq(Range_X[1],Range_X[2],length=ceiling(sqrt(Ngrid_sim)))
    Y = seq(Range_Y[1],Range_Y[2],length=ceiling(sqrt(Ngrid_sim)))
    # Expected recruitment
    image( x=X, y=Y, z=matrix(Omega_s[-c(1:(n_s+n_samp_per_year))]+log_mu_alpha, ncol=ceiling(sqrt(Ngrid_sim)), nrow=ceiling(sqrt(Ngrid_sim))), zlim=Zlim, col=Col(10), xaxt="n", yaxt="n" )
    mtext(side=3, line=0.5, outer=FALSE, text="Simulated")
    mtext(side=2, line=2, outer=FALSE, text=expression(italic(beta)+bold(A)))
    axis(2)
    image( x=X, y=Y, z=matrix(Report$Omega[loc_grid[,'Station_j']]+opt$par['beta'], ncol=ceiling(sqrt(Ngrid_sim)), nrow=ceiling(sqrt(Ngrid_sim))), zlim=Zlim, col=Col(10), xaxt="n", yaxt="n" )
    mtext(side=3, line=0.5, outer=FALSE, text="Estimated")
    # Realized recruitment
    for(t in 1:length(Years2Show)){
      image( x=X, y=Y, z=matrix(log(R_st[-c(1:(n_s+n_samp_per_year)),Years2Show[t]]), ncol=ceiling(sqrt(Ngrid_sim)), nrow=ceiling(sqrt(Ngrid_sim))), zlim=Zlim, col=Col(10), xaxt="n", yaxt="n" )
      if(t==length(Years2Show)) axis(1)
      axis(2)
      #if(t==1) mtext(side=3, line=0.5, outer=FALSE, text="True")
      mtext(side=2, line=2, outer=FALSE, text=switch(t,expression(bold(R[1])),expression(bold(R[2])),expression(bold(R[3]))) )
      image( x=X, y=Y, z=matrix(log(Report$R_it[cbind(loc_grid[,'Station_j'],Years2Show[t])]), ncol=ceiling(sqrt(Ngrid_sim)), nrow=ceiling(sqrt(Ngrid_sim))), zlim=Zlim, col=Col(10), xaxt="n", yaxt="n" )
      if(t==length(Years2Show)) axis(1)
      #if(t==1) mtext(side=3, line=0.5, outer=FALSE, text="Est")
    }
    mtext(side=1, outer=TRUE, line=2, text="Eastings (km)")
    mtext(side=4, outer=TRUE, line=0.5, text="Northings (km)")      
  dev.off()
  
  # Population trends
  png( file=paste(RunFile,"Est_S_and_N.png",sep=""), width=6, height=6, res=200, units="in")
    par(mfrow=c(2,2), mar=c(3,3,2,0))
    # Abundance
    Glm = glm( I_j ~ 0 + factor(Year_j) + factor(Station_j), data=DF, family="poisson")
    Mat = cbind( "Mean"=tapply(DF[,'I_j'],INDEX=DF[,'Year_j'],FUN=mean), "Glm"=exp(Glm$coef[1:length(unique(DF[,'Year_j']))]), "Year"=unique(DF[,'Year_j']))
    Mat[,c('Mean','Glm')] = Mat[,c('Mean','Glm')] * outer( rep(1,nrow(Mat)), mean(Report$N_t_hat)/c(mean(Mat[,"Mean"]),mean(Mat[,"Glm"])) ) 
    matplot( y=Report$N_t_hat, x=Year_Range[1]:Year_Range[2], type="l", col="black", lty="solid", ylim=c(0,max(Report$N_t_hat,Mat[,c('Mean','Glm')])), main="Abundance")
    matplot( y=Mat[,c("Mean",'Glm')], x=(Year_Range[1]:Year_Range[2])[Mat[,"Year"]], col=c("red","blue"), type="p", pch=21, add=TRUE)
    # Biomass
    Glm = glm( I(I_j*W_j) ~ 0 + factor(Year_j) + factor(Station_j), data=DF, family="gaussian")
    Mat = cbind( "Mean"=tapply(DF[,'W_j']*DF[,'I_j'],INDEX=DF[,'Year_j'],FUN=mean,na.rm=TRUE), "Glm"=Glm$coef[1:length(unique(DF[,'Year_j']))], "Year"=unique(DF[,'Year_j']))
    Mat[,c('Mean','Glm')] = Mat[,c('Mean','Glm')] * outer( rep(1,nrow(Mat)), mean(Report$S_t_hat)/c(mean(Mat[,"Mean"]),mean(Mat[,"Glm"])) ) 
    matplot( y=Report$S_t_hat, x=Year_Range[1]:Year_Range[2], type="l", col="black", lty="solid", ylim=c(0,max(Report$S_t_hat,Mat[,c('Mean','Glm')])), main="Spawning biomass")
    matplot( y=Mat[,c("Mean",'Glm')], x=(Year_Range[1]:Year_Range[2])[Mat[,"Year"]], col=c("red","blue"), type="p", pch=21, add=TRUE)
    # Average weight
    Glm = glm( I(ifelse(W_j==Inf,NA,W_j)) ~ 0 + factor(Year_j) + factor(Station_j), data=DF, family="gaussian", na.action="na.omit")
    Mat = cbind( "Mean"=tapply(DF[which(DF[,'I_j']>0),'W_j'],INDEX=DF[which(DF[,'I_j']>0),'Year_j'],FUN=mean,na.rm=TRUE), "Glm"=Glm$coef[1:length(unique(DF[,'Year_j']))], "Year"=unique(DF[,'Year_j']))
    Mat[,c('Mean','Glm')] = Mat[,c('Mean','Glm')] * outer( rep(1,nrow(Mat)), mean(Report$S_t_hat/Report$N_t_hat)/c(mean(Mat[,"Mean"]),mean(Mat[,"Glm"])) ) 
    matplot( y=Report$S_t_hat/Report$N_t_hat, x=Year_Range[1]:Year_Range[2], type="l", col="black", lty="solid", ylim=c(0,max(Report$S_t_hat/Report$N_t_hat,Mat[,c('Mean','Glm')])), main="Average weight")
    matplot( y=Mat[,c("Mean",'Glm')], x=(Year_Range[1]:Year_Range[2])[Mat[,"Year"]], col=c("red","blue"), type="p", pch=21, add=TRUE)
    # Recruitment
    matplot( y=Report$R_t_hat, x=Year_Range[1]:Year_Range[2], type="l", col="black", lty="solid", ylim=c(0,max(Report$R_t_hat)), main="Recruits")
  dev.off()

  # Fishing mortality
  png( file=paste(RunFile,"Est_Frate.png",sep=""), width=4, height=4, res=200, units="in")
    par(mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs="i")
    matplot( y=cbind(Report$F_t,Report$Exploit_Frac), x=Year_Range[1]:Year_Range[2], type="l", col="black", lty=c("solid","dotted"), ylim=c(0,max(cbind(Report$F_t,Report$Exploit_Frac))), xlab="", ylab="")
  dev.off()
  
  # Catch history
  png( file=paste(RunFile,"Pred_v_Obs_Catch.png",sep=""), width=4, height=4, res=200, units="in")
    matplot( y=cbind(C_t,Report$C_t_hat), x=Year_Range[1]:Year_Range[2], col=c("black","red"), type=c("l","p"), pch=21, ylim=c(0,max(cbind(C_t,Report$C_t_hat))))
  dev.off()

  # Weight -- pred vs. obs
  png( file=paste(RunFile,"Pred_v_Obs_MeanWeight.png",sep=""), width=4, height=4, res=200, units="in")
    Mat = cbind( "Obs"=DF[which(DF[,'I_j']>0),'W_j'], "Pred"=Report$W_it[as.matrix(DF[which(DF[,'I_j']>0),c('Station_j','Year_j')])] )
    # Cbind( Report$W_it[as.matrix(DF[which(DF[,'I_j']>0),c('Station_j','Year_j')])], (Report$S_it/Report$N_it)[as.matrix(DF[which(DF[,'I_j']>0),c('Station_j','Year_j')])]
    plot(x=Mat[,"Pred"], y=Mat[,"Obs"], xlim=c(0,max(Mat,na.rm=TRUE)), ylim=c(0,max(Mat,na.rm=TRUE)), col=rgb(0,0,0,alpha=0.02), , ylab="Obs", xlab="Pred"  ) 
    abline(a=0, b=1, col="red")                                                         # site 42 is wrong somehow
  dev.off()
  
  # Catch rates -- pred vs. obs
  png( file=paste(RunFile,"Pred_v_Obs_CatchRates.png",sep=""), width=4, height=4, res=200, units="in")
    Mat = cbind( "Obs"=DF[,'I_j'], "Pred"=Report$I_j_hat )
    plot(x=sqrt(Mat[,"Pred"]), y=sqrt(Mat[,"Obs"]), xlim=c(0,max(sqrt(Mat),na.rm=TRUE)), ylim=c(0,max(sqrt(Mat),na.rm=TRUE)), xaxt="n", yaxt="n", col=rgb(0,0,0,alpha=0.02), ylab="Obs", xlab="Pred"  ) 
    axis(1, at=axTicks(1), labels=axTicks(1)^2)
    axis(2, at=axTicks(2), labels=axTicks(2)^2)
    abline(a=0, b=1, col="red")
  dev.off()
              
  ##### Fields
  # Productivity
  png( file=paste(RunFile,"Fields_Omega_LL.png",sep=""), width=3, height=6, res=200, units="in")
    par(mfrow=c(2,1), mar=c(0,0,0,0))
    plot(1, type="n", ylim=range(loc_stations[,'long']), xlim=range(loc_stations[,'lat']))
    Bin = Bin_Quantile( Report$Omega[1:n_s], Nregions=5 )  
    points(y=loc_stations[,'lat'], x=loc_stations[,'long'], col=Col(Bin$Nregions)[Bin$Region], cex=3, pch=20)
    Legend( Bin, Col )
  dev.off()
  # Dynamical states
  Ncol = ceiling(sqrt(n_t+1)); Nrow = ceiling( (n_t+1)/Ncol )
  for(FigI in 1:5){
    if(FigI==1) Mat = (Report$Epsilon)[1:n_s,]
    if(FigI==2) Mat = log(Report$S_it)[1:n_s,]
    if(FigI==3) Mat = log(Report$W_it)[1:n_s,]
    if(FigI==4) Mat = log(Report$N_it)[1:n_s,]
    if(FigI==5) Mat = log(Report$R_it)[1:n_s,]
    png( file=paste(RunFile,"Fields_",c("Epsilon","S","W","N","R")[FigI],"_LL.png",sep=""), width=4*Ncol, height=2*Nrow, res=200, units="in")
      par(mfrow=c(Nrow,Ncol), mar=c(0,0,0,0), oma=c(2,2,0,0) )
      for(i in 1:n_t){
        #map("worldHires", ylim=y.lim,xlim=x.lim,col="grey90",fill=T, main="", mar=c(0,0,2.5,0))
        plot(1, type="n", ylim=range(loc_stations[,'long']), xlim=range(loc_stations[,'lat']))
        title( (Year_Range[1]:Year_Range[2])[i] )
        #points(y=loc_stations[,'lat'], x=loc_stations[,'long'], col=Col(n=10)[ceiling(f(Mat[1:n_s,])[,i]*9)+1], cex=3, pch=20)
        Bin = Bin_Quantile( Mat, Nregions=5 )  
        points(y=loc_stations[,'lat'], x=loc_stations[,'long'], col=Col(Bin$Nregions)[Bin$Region[,i]], cex=3, pch=20)
      	if( (i-1)%%Ncol == 0) axis(2)
      	if( (i-1)/Ncol >= (Nrow-1)) axis(1,las=2)
      	box(bty="o",lwd=2)
      }
      Legend( Bin, Col )
    dev.off()
  }

  # Spatial residuals (red is positive bias; blue is negative bias)
  if(ModelSet[ModelI]!="Index"){
    attach( TmbList )
    Ncol = ceiling(sqrt(length(unique(DF[,'Year_j'])))); Nrow = ceiling( length(unique(DF[,'Year_j']))/Ncol )
    for(FigI in 1:2){
      if(FigI==1){
        Mat1 = tapply( (Report[["W_j_hat"]]-Data[["W_j"]])^2/(Report[["W_j_hat"]]*(Data$CV_w+Report[["extraCV_w"]]))^2, INDEX=list(Data[['Station_j']], Data[['Year_j']]), FUN=mean, na.rm=TRUE)
        Mat2 = tapply( (Report[["W_j_hat"]]-Data[["W_j"]]), INDEX=list(Data[['Station_j']], Data[['Year_j']]), FUN=mean, na.rm=TRUE)
      }
      if(FigI==2){
        Mat1 = tapply( (Report[["I_j_hat"]]-Data[["I_j"]])^2/Report[["I_j_var"]], INDEX=list(Data[['Station_j']], Data[['Year_j']]), FUN=mean)
        Mat2 = tapply( (Report[["I_j_hat"]]-Data[["I_j"]]), INDEX=list(Data[['Station_j']], Data[['Year_j']]), FUN=mean)
      }
      # lat/lon Pearson residuals
      png( file=paste(RunFile,"PearsonResid_",c("W","I")[FigI],"_LL_stations.png",sep=""), width=4*Ncol, height=2*Nrow, res=200, units="in")
        par(mfrow=c(Nrow,Ncol), mar=c(0,0,0,0), oma=c(2,2,0,0) )
        for(i in 1:length(unique(DF[,'Year_j']))){
          #map("worldHires", ylim=y.lim,xlim=x.lim,col="grey90",fill=T, main="", mar=c(0,0,2.5,0))
          plot(1, type="n", ylim=range(loc_stations[,'long']), xlim=range(loc_stations[,'lat']))
          title( unique(DF[,'Year_j'])[i] )
          points(y=loc_stations[,'lat'], x=loc_stations[,'long'], col=c("blue","red")[ifelse(Mat2[,i]>0,2,1)], cex=sqrt(Mat1[,i]), pch=21)
          points(y=loc_stations[,'lat'], x=loc_stations[,'long'], col="black", cex=sqrt(0.5), pch=20)
        	if( (i-1)%%Ncol == 0) axis(2)
        	if( (i-1)/Ncol >= (Nrow-1)) axis(1,las=2)
        	box(bty="o",lwd=2)
        }
      dev.off()
      # individual Pearson residuals
      png( file=paste(RunFile,"PearsonResid_",c("W","I")[FigI],".png",sep=""), width=3*Ncol, height=3*Nrow, res=200, units="in")
        par(mfrow=c(Nrow,Ncol), mar=c(3,3,2,0) )
        for(i in 1:length(unique(DF[,'Year_j']))){
          Which = which(Data[["Year_j"]]==unique(Data[["Year_j"]])[i])
          if(FigI==1){
            Vec1 = ( (Report[["W_j_hat"]]-Data[["W_j"]])^2/(Report[["W_j_hat"]]*(Data$CV_w+Report[["extraCV_w"]]))^2 )
            Vec2 = ( (Report[["W_j_hat"]]-Data[["W_j"]]) )
          }
          if(FigI==2){
            Vec1 = ( (Report[["I_j_hat"]]-Data[["I_j"]])^2/Report[["I_j_var"]] )
            Vec2 = ( (Report[["I_j_hat"]]-Data[["I_j"]]) )
          }
          plot(y=ifelse(Vec2[Which]>0,1,-1)*Vec1[Which], x=Data[["Station_j"]][Which], pch=21, main=unique(Data[["Year_j"]])[i], ylim=RangeFn(ifelse(Vec2>0,1,-1)*Vec1))
        }
      dev.off()
      # individual Pearson residuals on map
      png( file=paste(RunFile,"PearsonResid_",c("W","I")[FigI],"_LL.png",sep=""), width=4*Ncol, height=2*Nrow, res=200, units="in")
        par(mfrow=c(Nrow,Ncol), mar=c(3,3,2,0) )
        for(i in 1:length(unique(DF[,'Year_j']))){
          Which = which(Data[["Year_j"]]==unique(Data[["Year_j"]])[i])
          if(FigI==1){
            Vec1 = ( (Report[["W_j_hat"]]-Data[["W_j"]])^2/(Report[["W_j_hat"]]*(Data$CV_w+Report[["extraCV_w"]]))^2 )
            Vec2 = ( (Report[["W_j_hat"]]-Data[["W_j"]]) )
          }
          if(FigI==2){
            Vec1 = ( (Report[["I_j_hat"]]-Data[["I_j"]])^2/Report[["I_j_var"]] )
            Vec2 = ( (Report[["I_j_hat"]]-Data[["I_j"]]) )
          }
          #map("worldHires", ylim=y.lim,xlim=x.lim,col="grey90",fill=T, main="", mar=c(0,0,2.5,0))
          plot(1, type="n", ylim=range(loc_stations[,'long']), xlim=range(loc_stations[,'lat']))
          title( unique(DF[,'Year_j'])[i] )
          points(y=DF[Which,'lat_j'], x=DF[Which,'long_j'], col=c("blue","red")[ifelse(Vec2[Which]>0,2,1)], cex=sqrt(Vec1[Which]), pch=21)
          points(y=DF[Which,'lat_j'], x=DF[Which,'long_j'],, col="black", cex=sqrt(0.5), pch=20)
        	if( (i-1)%%Ncol == 0) axis(2)
        	if( (i-1)/Ncol >= (Nrow-1)) axis(1,las=2)
        	box(bty="o",lwd=2)
          #plot(y=ifelse(Vec2[Which]>0,1,-1)*Vec1[Which], x=Data[["Station_j"]][Which], pch=21, main=unique(Data[["Year_j"]])[i], ylim=RangeFn(ifelse(Vec2>0,1,-1)*Vec1))
        }
      dev.off()
    } # End FigI loop
    detach( TmbList )
  } # End if(ModelI!="Index")

  # Detach data
  detach(DataList)
       