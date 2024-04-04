#########################################################################################################
# Comparing TESS models for simulated trees
#########################################################################################################

#Load necessary packages
x <- c("ape", "TESS", "diversitree", "geiger", "TreeSim", "paleotree", "dplyr")
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

#Directly simulate mass extinction with TreeSim (10% of species survive extinction at mass_ext_time)
#mass_ext_phy <- NULL
#mass_ext_time <- 66 #66 mya
#root_age <- 0
#while(is.null(mass_ext_phy) || isFALSE(between(root_age, 130, 140))){
#  mass_ext_phy <- sim.rateshift.taxa(n=70000, numbsim=1, lambda=c(1,2), mu=c(0.6, 0.8), frac=c(1,0.1), times=c(0, mass_ext_time))
#  mass_ext_extant <- drop.extinct(mass_ext_phy[[1]])
#  node_ages <- suppressMessages(dateNodes(mass_ext_extant))
#  node_ages <- sort(node_ages, decreasing = T)
#  root_age <- node_ages[1]
#}
#mass_ext_extant <- drop.extinct(mass_ext_phy[[1]])

#Check that the mass extinction occurs within the age of the tree
#node_ages <- dateNodes(mass_ext_extant)
#node_ages <- sort(node_ages, decreasing = T)
#node_ages[1]

#ALTERNATIVELY:
#mass_ext_phy <- NULL ; mass_ext_time <- 66 ; root_age <- 0 ; i=0
#while(is.null(mass_ext_phy) || isFALSE(between(root_age, 130, 140))){
#  lambda1 <- round(runif(n=1, min=0.1, max=1.5), 1)
#  lambda2 <- lambda1 * round(runif(n=1, min=1.5, max=2.5), 1)
#  mu1 <- lambda1 / round(runif(n=1, min=2, max=4), 1)
#  mu2 <- mu1 * round(runif(n=1, min=1.2, max=2), 1)
#  mass_ext_phy <- sim.rateshift.taxa(n=70000, numbsim=1, lambda=c(lambda1,lambda2), mu=c(mu1, mu2), frac=c(1,0.1), times=c(0, mass_ext_time))
#  mass_ext_extant <- drop.extinct(mass_ext_phy[[1]])
#  node_ages <- suppressMessages(dateNodes(mass_ext_extant))
#  node_ages <- sort(node_ages, decreasing = T)
#  root_age <- node_ages[1]
#  i = i+1
#  print(paste0("This while loop has so far run ", i, " iterations without success."))
#}

mass_ext_phy <- NULL
fraction <- 0.3
while(is.null(mass_ext_phy)){
  mass_ext_phy <- sim.rateshift.taxa(n=70000, numbsim=5, lambda=c(0.2,0.3), mu=c(0.1, 0.25), frac=c(1,fraction), times=c(0, 66))
}
for(i in 1:length(mass_ext_phy)){
  node_ages <- suppressMessages(dateNodes(mass_ext_phy[[i]]))
  node_ages <- sort(node_ages, decreasing = T)
  node_ages[1]
  print(paste0("Tree ", i, " has ", length(mass_ext_phy[[i]]$tip.label), " tips and a root age of ", node_ages[1], "."))
}

mass_ext_tree <- mass_ext_phy[[3]]
mass_ext_extant <- drop.extinct(mass_ext_tree)

#########################################################################################################
# Running the model comparison function on our simulated trees
#########################################################################################################
samplingfrac=0.22
extinct_removed_models <- model_comparison(sampfrac=samplingfrac, phy=extinct_removed)
mass_ext_extant_models <- model_comparison(sampfrac=samplingfrac, phy=mass_ext_extant)
#Constant BD strongly favored in both

#Now let's make some figures, using the mass extinction simulated phylogeny
tess.analysis(mass_ext_extant, empiricalHyperPriors = TRUE, samplingProbability = samplingfrac,
              estimateNumberMassExtinctions = TRUE, MAX_ITERATIONS = 10000, dir = "comet_mass_ext")
output <- tess.process.output("comet_mass_ext", numExpectedRateChanges = 1, numExpectedMassExtinctions = 1)
layout.mat <- matrix(1:6, nrow=3, ncol=2, byrow=TRUE)
layout(layout.mat)
tess.plot.output(output, fig.types = c("speciation rates", "speciation shift times", "extinction rates", "extinction shift times", "mass extinction Bayes factors", "mass extinction times"), las=2)

#Another attempt, setting hyperpriors manually rather than estimating them empirically:
#speciationPriorMu <- 0.2
#speciationPriorSigma <- 0.5
#extinctionPriorMu <- 0.1
#extinctionPriorSigma <- 0.5
#speciationRatePriorMean <- log((speciationPriorMu^2) / sqrt(speciationPriorSigma^2 + speciationPriorMu^2))
#speciationRatePriorStDev <- sqrt(log(1+speciationPriorSigma^2 / (speciationPriorMu^2)))
#extinctionRatePriorMean <- log((extinctionPriorMu^2) / sqrt(extinctionPriorSigma^2 + extinctionPriorMu^2))
#extinctionRatePriorStDev <- sqrt(log(1 + extinctionPriorSigma^2 / (extinctionPriorMu^2)))
#expectedSurvivalProbability <- 0.1
#pMassExtinctionPriorShape2 <- 100
#pMassExtinctionPriorShape1 <- -(pMassExtinctionPriorShape2 * expectedSurvivalProbability) / (expectedSurvivalProbability - 1)
#tess.analysis(mass_ext_extant, empiricalHyperPriors = FALSE, initialSpeciationRate = speciationPriorMu, 
#              speciationRatePriorMean = speciationRatePriorMean, speciationRatePriorStDev = speciationRatePriorStDev, 
#              initialExtinctionRate = extinctionPriorMu, extinctionRatePriorMean = extinctionRatePriorMean, 
#              extinctionRatePriorStDev = extinctionRatePriorStDev, samplingProbability = samplingfrac, 
#              numExpectedRateChanges = 1, numExpectedMassExtinctions = 1,
#              pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1, pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2, 
#              MAX_ITERATIONS = 10000, dir = "comet_mass_ext_two")
#output <- tess.process.output("comet_mass_ext_two", numExpectedRateChanges = 2, numExpectedMassExtinctions = 2)
#layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
#layout(layout.mat)
#tess.plot.output(output, fig.types = c("speciation rates", "speciation shift times", "extinction rates", "extinction shift times", "mass extinction Bayes factors", "mass extinction times"), las=2)

#REVIEWER SUGGESTION: Grand Coupure dampening?
two_mass_exts <- NULL
fractions <- c(0.1, 0.8)
while(is.null(two_mass_exts)){
 #33.5 mya is Eocene-Oligocene boundary
 two_mass_exts <- sim.rateshift.taxa(n=70000, numbsim=5, lambda=c(0.15,0.2,0.25), mu=c(0.13, 0.18, 0.22), frac=c(1,fractions), times=c(0, 33.5, 66))
}
for(i in 1:length(two_mass_exts)){
 node_ages <- suppressMessages(dateNodes(two_mass_exts[[i]]))
 node_ages <- sort(node_ages, decreasing = T)
 node_ages[1]
 print(paste0("Tree ", i, " has ", length(two_mass_exts[[i]]$tip.label), " tips and a root age of ", node_ages[1], "."))
}

two_mass_exts_tree <- two_mass_exts[[2]]
two_mass_exts_extant <- drop.extinct(two_mass_exts_tree)

tess.analysis(two_mass_exts_extant, empiricalHyperPriors = TRUE, samplingProbability = samplingfrac,
              estimateNumberMassExtinctions = TRUE, MAX_ITERATIONS = 10000, dir = "comet_mass_ext")
output <- tess.process.output("comet_mass_ext", numExpectedRateChanges = 2, numExpectedMassExtinctions = 2)
layout.mat <- matrix(1:6, nrow=3, ncol=2, byrow=TRUE)
layout(layout.mat)
tess.plot.output(output, fig.types = c("speciation rates", "speciation shift times", "extinction rates", "extinction shift times", "mass extinction Bayes factors", "mass extinction times"), las=2)



