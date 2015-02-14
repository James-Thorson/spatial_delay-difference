\name{MakeTmbList_Fn}
\alias{MakeTmbList_Fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a tagged list of inputs for TMB
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MakeTmbList_Fn(Version, Model, Fix_Q, ErrorModel_CatchRates, ErrorModel, Smooth_F, n_j, n_i, n_s, n_t, DF_input, C_t, mesh_stations, spde_stations, Area_i, alpha_g, ro, w_k, M, k, CV_c, CV_w, q_I)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Version}{
%%     ~~Describe \code{Version} here~~
}
  \item{Model}{
%%     ~~Describe \code{Model} here~~
}
  \item{Fix_Q}{
%%     ~~Describe \code{Fix_Q} here~~
}
  \item{ErrorModel_CatchRates}{
%%     ~~Describe \code{ErrorModel_CatchRates} here~~
}
  \item{ErrorModel}{
%%     ~~Describe \code{ErrorModel} here~~
}
  \item{Smooth_F}{
%%     ~~Describe \code{Smooth_F} here~~
}
  \item{n_j}{
%%     ~~Describe \code{n_j} here~~
}
  \item{n_i}{
%%     ~~Describe \code{n_i} here~~
}
  \item{n_s}{
%%     ~~Describe \code{n_s} here~~
}
  \item{n_t}{
%%     ~~Describe \code{n_t} here~~
}
  \item{DF_input}{
%%     ~~Describe \code{DF_input} here~~
}
  \item{C_t}{
%%     ~~Describe \code{C_t} here~~
}
  \item{mesh_stations}{
%%     ~~Describe \code{mesh_stations} here~~
}
  \item{spde_stations}{
%%     ~~Describe \code{spde_stations} here~~
}
  \item{Area_i}{
%%     ~~Describe \code{Area_i} here~~
}
  \item{alpha_g}{
%%     ~~Describe \code{alpha_g} here~~
}
  \item{ro}{
%%     ~~Describe \code{ro} here~~
}
  \item{w_k}{
%%     ~~Describe \code{w_k} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{CV_c}{
%%     ~~Describe \code{CV_c} here~~
}
  \item{CV_w}{
%%     ~~Describe \code{CV_w} here~~
}
  \item{q_I}{
%%     ~~Describe \code{q_I} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Version, Model, Fix_Q, ErrorModel_CatchRates, ErrorModel, 
    Smooth_F, n_j, n_i, n_s, n_t, DF_input, C_t, mesh_stations, 
    spde_stations, Area_i, alpha_g, ro, w_k, M, k, CV_c, CV_w, 
    q_I) 
{
    F_t = rep(0.1, n_t)
    if (Version == "delay_difference_v4d") {
        Error_Model = ErrorModel_CatchRates
        Data = list(Error_Model = Error_Model, Smooth_F = Smooth_F, 
            n_j = n_j, n_i = n_i, n_s = n_s, n_t = n_t, I_j = DF_input[, 
                "I_j"], W_j = DF_input[, "W_j"], AreaSwept_j = DF_input[, 
                "AreaSwept_j"], Station_j = DF_input[, "Station_j"] - 
                1, Year_j = DF_input[, "Year_j"] - 1, C_t = C_t, 
            meshidxloc = mesh_stations$idx$loc - 1, G0 = spde_stations$param.inla$M0, 
            G1 = spde_stations$param.inla$M1, G2 = spde_stations$param.inla$M2, 
            Area_i = Area_i, alpha_g = alpha_g, ro = ro, w_k = w_k, 
            M = M, k = k, CV_c = CV_c, CV_w = CV_w)
        Parameters = list(log_F_sd = log(1), log_F_equil = log(F_t[1]), 
            log_F_t_input = log(F_t), log_q_I = log(q_I), beta = c(0), 
            log_tau_E = 0, log_tau_O = 0, log_kappa = 0, ln_VarInfl = c(0, 
                0), Epsilon_input = matrix(0, spde_stations$n.spde, 
                n_t), Omega_input = rep(0, spde_stations$n.spde))
    }
    if (Version == "delay_difference_v6c") {
        Data = list(ModelType = ifelse(Model == "Spatial", 1, 
            2), ErrorModel_CatchRates = ErrorModel_CatchRates, 
            ErrorModel_MeanWeight = ErrorModel_MeanWeight, Smooth_F = Smooth_F, 
            n_j = n_j, n_i = n_i, n_s = n_s, n_t = n_t, I_j = DF_input[, 
                "I_j"], W_j = DF_input[, "W_j"], AreaSwept_j = DF_input[, 
                "AreaSwept_j"], Station_j = DF_input[, "Station_j"] - 
                1, Year_j = DF_input[, "Year_j"] - 1, C_t = C_t, 
            meshidxloc = mesh_stations$idx$loc - 1, G0 = spde_stations$param.inla$M0, 
            G1 = spde_stations$param.inla$M1, G2 = spde_stations$param.inla$M2, 
            Area_i = Area_i, alpha_g = alpha_g, ro = ro, w_k = w_k, 
            M = M, k = k, CV_c = CV_c, CV_w = CV_w)
        Parameters = list(log_F_sd = log(1), log_F_equil = log(F_equil), 
            log_F_t_input = log(F_t), log_q_I = log(q_I), beta = c(0), 
            log_tau_E = 0, log_tau_O = 0, log_kappa = 0, ln_VarInfl = c(0, 
                0), log_extraCV_w = log(0.05), log_tau_N = log(1), 
            Epsilon_input = matrix(0, spde_stations$n.spde, n_t), 
            Omega_input = rep(0, spde_stations$n.spde), Nu_input = rep(0, 
                n_t))
    }
    if (Version == "delay_difference_v8d") {
        Data = list(ModelType = ifelse(Model == "Spatial", 1, 
            2), ErrorModel_CatchRates = ErrorModel_CatchRates, 
            ErrorModel_MeanWeight = ErrorModel_MeanWeight, Smooth_F = Smooth_F, 
            n_j = n_j, n_i = n_i, n_s = n_s, n_t = n_t, I_j = DF_input[, 
                "I_j"], W_j = DF_input[, "W_j"], AreaSwept_j = DF_input[, 
                "AreaSwept_j"], Station_j = DF_input[, "Station_j"] - 
                1, Year_j = DF_input[, "Year_j"] - 1, C_t = C_t, 
            IndexMat = IndexMat, meshidxloc = mesh_stations$idx$loc - 
                1, G0 = spde_stations$param.inla$M0, G1 = spde_stations$param.inla$M1, 
            G2 = spde_stations$param.inla$M2, Area_i = Area_i, 
            alpha_g = alpha_g, ro = ro, w_k = w_k, M = M, k = k, 
            CV_c = CV_c, CV_w = CV_w)
        Parameters = list(log_F_sd = log(1), log_F_t_input = log(c(F_equil, 
            F_t)), log_q_I = log(q_I), beta = log_mu_alpha, log_tau_E = 0, 
            log_tau_O = 0, log_kappa = 0, ln_VarInfl = c(0, 0), 
            log_extraCV_w = log(0.05), log_tau_N = log(1), log_extraCV_Index = rep(log(0.1), 
                2), Epsilon_input = matrix(0, spde_stations$n.spde, 
                n_t), Omega_input = rep(0, spde_stations$n.spde), 
            Nu_input = rep(0, n_t))
    }
    if (Model == "Spatial") {
        if (Smooth_F == 0) 
            Random = c("Epsilon_input", "Omega_input")
        if (Smooth_F != 0) 
            Random = c("Epsilon_input", "Omega_input", "log_F_equil")
    }
    if (Model == "Nonspatial") {
        if (Smooth_F == 0) 
            Random = c("Nu_input")
        if (Smooth_F != 0) 
            Random = c("Nu_input", "log_F_t_input")
    }
    if (Model == "Index") {
        if (Smooth_F == 0) 
            Random = c("Nu_input")
        if (Smooth_F != 0) 
            Random = c("Nu_input", "log_F_t_input")
    }
    Map = list()
    if (Smooth_F == 0) 
        Map[["log_F_sd"]] = factor(NA)
    if (ErrorModel_CatchRates == 0) 
        Map[["ln_VarInfl"]] = factor(c(NA, NA))
    if (ErrorModel_MeanWeight == 0) 
        Map[["log_extraCV_w"]] = factor(NA)
    if (Model == "Spatial") {
        Map[["log_tau_N"]] = factor(NA)
        Map[["Nu_input"]] = factor(rep(NA, n_t))
        Map[["log_extraCV_Index"]] = factor(rep(NA, 2))
        Map[["log_F_t_input"]] = factor(c(1, 1:n_t))
    }
    if (Model == "Nonspatial") {
        Map[["log_tau_E"]] = factor(NA)
        Map[["log_tau_O"]] = factor(NA)
        Map[["Epsilon_input"]] = factor(rep(NA, length(Parameters$Epsilon_input)))
        Map[["Omega_input"]] = factor(rep(NA, length(Parameters$Omega_input)))
        Map[["log_kappa"]] = factor(NA)
        Map[["log_extraCV_Index"]] = factor(rep(NA, 2))
        Map[["log_F_t_input"]] = factor(c(1, 1:n_t))
    }
    if (Model == "Index") {
        Map[["log_tau_E"]] = factor(NA)
        Map[["log_tau_O"]] = factor(NA)
        Map[["Epsilon_input"]] = factor(rep(NA, length(Parameters$Epsilon_input)))
        Map[["Omega_input"]] = factor(rep(NA, length(Parameters$Omega_input)))
        Map[["log_kappa"]] = factor(NA)
        Map[["ln_VarInfl"]] = factor(c(NA, NA))
        Map[["log_extraCV_w"]] = factor(NA)
        Map[["log_F_t_input"]] = factor(c(1, 1:n_t))
    }
    if (Fix_Q == TRUE) 
        Map[["log_q_I"]] = factor(NA)
    Return = list(Parameters = Parameters, Data = Data, Map = Map, 
        Random = Random)
    return(Return)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
