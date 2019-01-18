# Fit Bayesian univariate with Stan
# Created  ~Nov-2018
# Shifted between folders/projects for principled expansion and tidied Jan-17-2019

# Set up ----------------------------------------------------------------------
rm(list=ls()) # Clears everything from environment

# Below: function to load necessary packages and install if necessary
init <- function(need) {
  ip <- .packages(all.available = T)
  if (any((need %in% ip) == F)) {
    install.packages(need[!(need %in% ip)])
  }
  ok <- sapply(1:length(need), function(p) require(need[[p]], 
                                                   character.only = T))
}

init(c("rstan","plyr","reshape2"))
options(mc.cores = parallel::detectCores())

# Load files ----------------------------------------------------

source("Preprocessing_data.R")

# Parameters and inits

parameters <- c("theta_signalitems", "theta_noiseitems",
                "grand_mu", "grand_sigma",
                "dprime","bias","sigma_ratio",
                "crit", "mu_crits", "sigma_crits",
                "Omega_subj", "sigma_alpha_subj","alpha_subj",
                "lp__")

initsplease <- function(numberPara = 4) {
  list(
    grand_mu=sort(rnorm(2, 0.5), decreasing = TRUE),  
    grand_sigma = runif(2, 1, 2),
    mu_crits = sort(rnorm(3)), sigma_crits = runif(2),
    crit_subj = t(sapply(1:numberParticipants, function(x) sort(rnorm(3)))),
    Lower_Omega_subj=diag(numberPara ), 
    sigma_alpha_subj = runif(numberPara ), 
    alpha_t_subj=matrix(rnorm(numberParticipants * numberPara , 0, 0.1), 
                        numberPara , numberParticipants)
  )
}



fit_uvsdt_participant <- stan(file = "Univariate/UVSDT_participant.stan",   
                              data=cidrs18Item_info, 
                              init=initsplease,
                              pars=parameters,
                              iter=4000, 
                              chains=6, 
                              thin=3,
                              warmup = 1000
)


# Quick results
print(fit_uvsdt_participant, pars = c("dprime","bias","sigma_ratio"))
print(fit_uvsdt_participant, pars = c("grand_mu", "grand_sigma", "lp__"), include = TRUE)

# MCMC diagnostics
summary(fit_uvsdt_participant)$summary[,"Rhat"]

stan_ac(fit_uvsdt_participant, pars = c("grand_mu", "grand_sigma", "lp__"))
stan_ac(fit_uvsdt_participant, pars = c("Omega_subj","sigma_alpha_subj", "alpha_subj"))

rstan::traceplot(fit_uvsdt_participant, pars = c("grand_mu", "grand_sigma", "lp__"))
rstan::traceplot(fit_uvsdt_participant, pars = c("Omega_subj", "sigma_alpha_subj","lp__"))
rstan::traceplot(fit_uvsdt_participant, pars = c("alpha_subj", "lp__"))


save(fit_uvsdt_participant, fit_uvsdt_participant, file = "fit_uvsdt_participant.Rda", compress = "xz")

# Bit more involved results
# highest posterior density interval
getCredI <- function(x, interval = 0.95) {
  ci <- coda::HPDinterval(coda::as.mcmc(as.vector(x)), prob = interval)
  c(mean = mean(x), lower = ci[1], upper = ci[2])
}


  
someResults <- function(Model){
  post <- rstan::extract(Model, pars = c("grand_mu", "grand_sigma",
                                                       "alpha_subj","Omega_subj",
                                                       "mu_crits","sigma_crits",
                                                       "crit"))
  mu_diff <- post$grand_mu[,1] - post$grand_mu[,2]
  poolsigma<- sqrt( ( post$grand_sigma[,1]^2 + post$grand_sigma[,2]^2)/2)
  dsuba <- getCredI(mu_diff/poolsigma)
  sdtbias <- (post$grand_mu[,1] + post$grand_mu[,2])/2 #look up c_adj
  cbias <- getCredI(sdtbias)
  participant_effects <- apply(post$alpha_subj,2,getCredI)

  return(list(dsuba,cbias,participant_effects))
}


someResults(fit_uvsdt_participant)
