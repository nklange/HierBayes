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

# highest posterior density interval, is this the best function to do this?
getCredI <- function(x, interval = 0.95) {
  ci <- coda::HPDinterval(coda::as.mcmc(as.vector(x)), prob = interval)
  c(mean = mean(x), lower = ci[1], upper = ci[2])
}


# Unequal Variance ----------------------------------------------
# Parameters and inits

parameters <- c("theta_signalitems", "theta_noiseitems",
                "grand_mu", "grand_sigma",
                "dprime","bias","sigma_ratio","dprime_std",
                "crit", "mu_crits", "sigma_crits",
                "Omega_subj", "sigma_alpha_subj","alpha_subj",
                "lp__")

initsplease <- function(numberPara = 4) {
  list(
    grand_mu=sort(rnorm(2, 0.5), decreasing = TRUE),  
    grand_sigma = runif(2, 1, 2),
    mu_crits = sort(rnorm(3)), sigma_crits = runif(2),
    crit_un = t(sapply(1:numberParticipants, function(x) sort(rnorm(3)))),
    Lower_Omega_subj=diag(numberPara ), 
    sigma_alpha_subj = runif(numberPara ), 
    alpha_t_subj=matrix(rnorm(numberParticipants * numberPara , 0, 0.1), 
                        numberPara , numberParticipants)
  )
}

#Data:
## cidrs18Item_info
## cidrs18Context_info

fit_uvsdt_participant <- stan(file = "Univariate/UVSDT_participant.stan",   
                              data=cidrs18Context_info, 
                              init=initsplease,
                              pars=parameters,
                              iter=4000, 
                              chains=4, 
                              thin=3,
                              warmup = 1000
)

save(fit_uvsdt_participant, file = "cidrs18Context_uvsdt_participant.Rda", compress = "xz")


# Quick results
print(fit_uvsdt_participant, pars = c("dprime","bias","sigma_ratio","dprime_std"))
print(fit_uvsdt_participant, pars = c("grand_mu", "grand_sigma", "lp__"), include = TRUE)

# MCMC diagnostics
summary(fit_uvsdt_participant)$summary[,"Rhat"]

stan_ac(fit_uvsdt_participant, pars = c("grand_mu", "grand_sigma", "lp__"))
stan_ac(fit_uvsdt_participant, pars = c("Omega_subj","sigma_alpha_subj", "alpha_subj"))

rstan::traceplot(fit_uvsdt_participant, pars = c("grand_mu", "grand_sigma", "lp__"))
rstan::traceplot(fit_uvsdt_participant, pars = c("Omega_subj", "sigma_alpha_subj","lp__"))
rstan::traceplot(fit_uvsdt_participant, pars = c("alpha_subj", "lp__"))


someResults <- function(Model,numberPara = 4){
  post <- rstan::extract(Model, pars = c("grand_mu", "grand_sigma",
                                         "alpha_subj","Omega_subj",
                                         "mu_crits","sigma_crits",
                                         "crit"))
  mu_diff <- post$grand_mu[,1] - post$grand_mu[,2]
  poolsigma<- sqrt( ( post$grand_sigma[,1]^2 + post$grand_sigma[,2]^2)/2)
  dsuba <- getCredI(mu_diff/poolsigma)
  sdtbias <- (post$grand_mu[,1] + post$grand_mu[,2])/2 #look up c_adj
  cbias <- getCredI(sdtbias)

  participant_effects <- list()
  for (i in 1:numberPara){
  
    participant_effects[[i]]<-apply(post$alpha_subj[,,i],2,getCredI)
  
  }

  return(list(dsuba,cbias,participant_effects))
}


someResults(fit_uvsdt_participant,4)



# Equal Variance -----------------------------------------------------

# Parameters and inits

parameters <- c("theta_signalitems", "theta_noiseitems",
                "grand_mu", "grand_sigma",
                "dprime","bias",
                "crit", "mu_crits", "sigma_crits",
                "Omega_subj", "sigma_alpha_subj","alpha_subj",
                "lp__")

initsplease <- function(numberPara = 3) {
  list(
    grand_mu=sort(rnorm(2, 0.5), decreasing = TRUE),  
    grand_sigma = runif(1, 1, 2),
    mu_crits = sort(rnorm(3)), sigma_crits = runif(2),
    crit_un = t(sapply(1:numberParticipants, function(x) sort(rnorm(3)))),
    Lower_Omega_subj=diag(numberPara ), 
    sigma_alpha_subj = runif(numberPara ), 
    alpha_t_subj=matrix(rnorm(numberParticipants * numberPara , 0, 0.1), 
                        numberPara , numberParticipants)
  )
}

#Data:
## cidrs18Item_info
## cidrs18Context_info

fit_evsdt_participant <- stan(file = "Univariate/EVSDT_participant.stan",   
                              data=cidrs18Context_info, 
                              init=initsplease,
                              pars=parameters,
                              iter=4000, 
                              chains=6, 
                              thin=3,
                              warmup = 1000
)

save(fit_evsdt_participant, file = "cidrs18Context_evsdt_participant.Rda", compress = "xz")


# Quick results
print(fit_evsdt_participant, pars = c("dprime","bias"))
print(fit_evsdt_participant, pars = c("grand_mu", "grand_sigma", "lp__"), include = TRUE)

# MCMC diagnostics
summary(fit_evsdt_participant)$summary[,"Rhat"]

stan_ac(fit_evsdt_participant, pars = c("grand_mu", "grand_sigma", "lp__"))
stan_ac(fit_evsdt_participant, pars = c("Omega_subj","sigma_alpha_subj", "alpha_subj"))

rstan::traceplot(fit_evsdt_participant, pars = c("grand_mu", "grand_sigma", "lp__"))
rstan::traceplot(fit_evsdt_participant, pars = c("Omega_subj", "sigma_alpha_subj","lp__"))
rstan::traceplot(fit_evsdt_participant, pars = c("alpha_subj", "lp__"))





# Bit more involved results 

someResults <- function(Model, numberPara = 3){
  post <- rstan::extract(Model, pars = c("grand_mu", "grand_sigma",
                                         "alpha_subj","Omega_subj",
                                         "mu_crits","sigma_crits",
                                         "crit"))
  mu_diff <- post$grand_mu[,1] - post$grand_mu[,2]
  poolsigma<- sqrt(post$grand_sigma^2)
  dsuba <- getCredI(mu_diff/poolsigma)
  sdtbias <- (post$grand_mu[,1] + post$grand_mu[,2])/2 #look up c_adj
  cbias <- getCredI(sdtbias)
  
  participant_effects <- list()
  for (i in 1:numberPara){
    
    participant_effects[[i]]<-apply(post$alpha_subj[,,i],2,getCredI)
    
  }
  
  return(list(dsuba,cbias,participant_effects))
}


someResults(fit_evsdt_participant)



# DP Equal Variance -----------------------------------------------------

# Parameters and inits

parameters <- c("theta_signalitems", "theta_noiseitems",
                "grand_mu", "grand_sigma","grand_R",
                "familiarity",
                "crit", "mu_crits", "sigma_crits",
                "Omega_subj", "sigma_alpha_subj","alpha_subj",
                "lp__")

initsplease <- function(numberPara = 5) {
  list(
    grand_mu=sort(rnorm(2, 0.5), decreasing = TRUE),  
    grand_sigma = runif(1, 1, 2),
    grand_R = runif(2, 0, 1),
    mu_crits = sort(rnorm(3)), sigma_crits = runif(2),
    crit_un = t(sapply(1:numberParticipants, function(x) sort(rnorm(3)))),
    Lower_Omega_subj=diag(numberPara ), 
    sigma_alpha_subj = runif(numberPara ), 
    alpha_t_subj=matrix(rnorm(numberParticipants * numberPara , 0, 0.1), 
                        numberPara , numberParticipants)
  )
}

#Data:
## cidrs18Item_info
## cidrs18Context_info

fit_dpsource2R_participant <- stan(file = "Univariate/DP_EVSDT_Source2R_participant.stan",   
                              data=cidrs18ContextS_info, 
                              init=initsplease,
                              pars=parameters,
                              iter=4000, 
                              chains=6, 
                              thin=3,
                              warmup = 1000
)

save(fit_dpsource2R_participant, file = "cidrs18Context_fit_dpsource2R_participant_participant.Rda", compress = "xz")


# Quick results
print(fit_dpsource2R_participant, pars = c("familiarity"))
print(fit_dpsource2R_participant, pars = c("grand_mu", "grand_sigma","grand_R", "lp__"), include = TRUE)

# MCMC diagnostics
summary(fit_dpsource2R_participant)$summary[,"Rhat"]

stan_ac(fit_dpsource2R_participant, pars = c("grand_mu", "grand_sigma", "lp__"))
stan_ac(fit_dpsource2R_participant, pars = c("Omega_subj","sigma_alpha_subj", "alpha_subj"))

rstan::traceplot(fit_dpsource2R_participant, pars = c("grand_mu", "grand_sigma", "lp__"))
rstan::traceplot(fit_dpsource2R_participant, pars = c("Omega_subj", "sigma_alpha_subj","lp__"))
rstan::traceplot(fit_dpsource2R_participant, pars = c("alpha_subj", "lp__"))


# Bit more involved results 

someResults <- function(Model, numberPara = 5){
  post <- rstan::extract(Model, pars = c("grand_mu", "grand_sigma",
                                         "alpha_subj","Omega_subj",
                                         "mu_crits","sigma_crits",
                                         "crit"))
  mu_diff <- post$grand_mu[,1] - post$grand_mu[,2]
  poolsigma<- sqrt(post$grand_sigma^2)
  dsuba <- getCredI(mu_diff/poolsigma)
  
  participant_effects <- list()
  for (i in 1:numberPara){
    
    participant_effects[[i]]<-apply(post$alpha_subj[,,i],2,getCredI)
    
  }
  
  return(list(dsuba,participant_effects))
}


someResults(fit_dpsource2R_participant)



