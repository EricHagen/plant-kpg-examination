#########################################################################################################
# Comparing TESS models for simulated trees
#########################################################################################################

#Load necessary packages
x <- c("ape", "TESS", "diversitree", "geiger", "TreeSim")
invisible(lapply(x, library, character.only=T))

#Model comparison function
model_comparison <- function(sampfrac, phy){
  
  #Set necessary values
  times <- as.numeric(branching.times(phy)) #Branching times
  survivalProbability <- 0.1 #Survival probability for mass extinction model
  rateChangeTime <- max(times)/2 #Times of rate shifts in episodic model
  
  #Create some prior distributions
  prior_delta <- function(x){ 
    dexp(x, rate=10, log=TRUE) 
  }
  prior_delta_before <- prior_tau_before <- prior_delta_after <- prior_tau_after <- prior_tau <- prior_delta
  prior_time <- function(x){
    dunif(x, min=max(times)/2, max=max(times), log=TRUE)
  }
  
  #Write some prior functions
  priorsConstBD <- c("diversification"=prior_delta, "turnover"=prior_tau)
  priorsEpisodicBD <- c("diversification_before"=prior_delta_before, "turnover_before"=prior_tau_before, "diversification_after"=prior_delta_after, "turnover_after"=prior_tau_after)
  priorsMassExtinctionBD <- c("diversification"=prior_delta, "turnover"=prior_tau, "mass-extinctiontime"=prior_time)
  
  #Model functions
  constBD <- function(params){
    l <- params[1] + params[2]
    m <- params[2]
    lnl <- tess.likelihood(times, l, m, samplingProbability=sampfrac, samplingStrategy="diversified")
    return(lnl)
  }
  episodicBD <- function(params){
    speciation <- c(params[1]+params[2], params[3]+params[4])
    extinction <- c(params[2], params[4])
    lnl <- tess.likelihood.rateshift(times, lambda=speciation, mu=extinction, rateChangeTimesLambda=rateChangeTime, rateChangeTimesMu=rateChangeTime, samplingProbability=sampfrac, log=TRUE)
    return(lnl)
  }
  likelihoodMassExtinctionBD <- function(params) { 
    speciation <- params[1] + params[2] 
    extinction <- params[2] 
    time <- params[3] 
    lnl <- tess.likelihood(times, lambda=speciation, mu=extinction, massExtinctionTimes=time, massExtinctionSurvivalProbabilities=survivalProbability, samplingProbability=1.0, log=TRUE) 
    return(lnl) 
  }
  
  #Model fitting
  print("Fitting Constant BD (1 of 3)...")
  mlconstBD <- tess.steppingStoneSampling(likelihoodFunction=constBD, priors=priorsConstBD, parameters=runif(2,0,1), logTransforms=c(TRUE,TRUE), iterations=1000, burnin=100, K=50)
  print("Fitting Episodic BD (2 of 3)...")
  mlepisodicBD <- tess.steppingStoneSampling(likelihoodFunction = episodicBD, priors = priorsEpisodicBD, parameters = runif(4,0,1), logTransforms = c(TRUE,TRUE,TRUE,TRUE), iterations = 1000, burnin = 100, K = 50)
  print("Fitting Mass Extinction BD (3 of 3)...")
  mlMassExtinctionBD <- tess.steppingStoneSampling(likelihoodFunction=likelihoodMassExtinctionBD, priors=priorsMassExtinctionBD, parameters=c(runif(2,0,1), max(times)*3/4), logTransforms=c(TRUE,TRUE,FALSE), iterations=1000, burnin=100, K=50)
  
  #Organize & return data
  candidateModels <- c("ConstBD"=mlconstBD, "EpisodicBD"=mlepisodicBD, "MassExtBD"=mlMassExtinctionBD)
  marginalLikelihoodGrid <- expand.grid(M0=names(candidateModels), M1=names(candidateModels))
  marginalLikelihoodGrid$BF <- 2 * (candidateModels[marginalLikelihoodGrid$M0] - candidateModels[marginalLikelihoodGrid$M1])
  marginalLikelihoodGrid <- marginalLikelihoodGrid[order(marginalLikelihoodGrid$BF, decreasing=TRUE),]
  return(marginalLikelihoodGrid)
  
}

#########################################################################################################
# Simulating trees
#########################################################################################################

#Removing extinct taxa from a simulated phylogeny
extinct_removed_orig <- NULL
while(is.null(extinct_removed_orig) || length(extinct_removed_orig$tip.label) < 100){
  extinct_removed_orig <- tree.bd(pars=c(1, 0.3), max.taxa=10000, include.extinct=TRUE)
}
extinct_removed <- drop.extinct(extinct_removed_orig)

#Directly simulate mass extinction with TreeSim (10% of species survive extinction 30 mya)
mass_ext_phy <- NULL
while(is.null(mass_ext_phy)){
  mass_ext_phy <- sim.rateshift.taxa(n=5000, numbsim=1, lambda=c(1,0.5), mu=c(0.3, 0.4), frac=c(1,0.1), times=c(0,30))
}
mass_ext_extant <- drop.extinct(mass_ext_phy[[1]])

#########################################################################################################
# Running the model comparison function on our simulated trees
#########################################################################################################
samplingfrac=0.5
extinct_removed_models <- model_comparison(sampfrac=samplingfrac, phy=extinct_removed)
mass_ext_extant_models <- model_comparison(sampfrac=samplingfrac, phy=mass_ext_extant)


