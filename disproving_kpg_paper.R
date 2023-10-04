#########################################################################################################
# Disproving Thompson & Ramirez-Barahona 2023
#########################################################################################################

#Load necessary packages
x <- c("ape", "TreeTools", "phytools", "TESS", "diversitree", "geiger", "TreeSim", "paleotree")
invisible(lapply(x, library, character.only=T))

#Load necessary data
#sb_tree <- read.tree("/Users/erichagen1/Desktop/University_of_Arkansas_Stuff/Publications/polyploidy_and_extinctionrisk/tree_files/raw_smithandbrown.tre.txt")
sb_tree <- read.tree("raw_smithandbrown.tre.txt")

#For now we will just use the mesangiosperms for our test
sb_tree <- Preorder(sb_tree)
node_n <- which(sb_tree$node.label == "Mesangiospermae") + Ntip(sb_tree)
mesangiosperms <- Subtree(sb_tree, node_n) #78,821 taxa
sb_tree <- drop.tip(sb_tree, which(!sb_tree$tip.label %in% mesangiosperms$tip.label))

#LTT with phytools
pdf("ltt.pdf")
mesangiosperm_ltt <- ltt(sb_tree, plot=TRUE)
dev.off()

#Since the LTT looks weird, and force.ultrametric isn't working, let's try trimming the tree down just a bit
sb_tree_trimmed <- timeSliceTree(ttree=sb_tree, sliceTime=0.01, plot=F)
is.ultrametric(sb_tree_trimmed) #TRUE!!!
pdf("ltt_trimmed.pdf")
mesangiosperm_ltt_trimmed <- ltt(sb_tree, plot=TRUE)
dev.off()

#########################################################################################################
#Model comparison with TESS (made consulting suppl. from Chaitanya et al. 2023)
#########################################################################################################

#What is sampling fraction of mesangiosperms in the Smith & Brown tree?
sampfrac = length(sb_tree_trimmed$tip.label) / 350000 #Paton et al. 2008

#Branching times
times <- as.numeric(branching.times(sb_tree_trimmed))

#Constant birth-death
constBD <- function(params){
  l <- params[1] + params[2]
  m <- params[2]
  lnl <- tess.likelihood(times, l, m, samplingProbability=sampfrac, samplingStrategy="diversified")
  return(lnl)
}

prior_delta <- function(x){ 
  dexp(x, rate=10, log=TRUE) 
}
prior_tau <- function(x){
  dexp(x, rate=10, log=TRUE) 
}
priorsConstBD <- c("diversification"=prior_delta, "turnover"=prior_tau)

mlconstBD <- tess.steppingStoneSampling(likelihoodFunction=constBD, 
                              priors=priorsConstBD, parameters=runif(2,0,1), logTransforms=c(TRUE,TRUE), 
                              iterations=1000, burnin=100, K=50)

#Episodic birth-death on S&B tree
prior_delta_before <- function(x){ 
  dexp(x,rate=10.0,log=TRUE)
}
prior_tau_before <- function(x){ 
  dexp(x,rate=10.0,log=TRUE)
}
prior_delta_after <- function(x){ 
  dexp(x,rate=10.0,log=TRUE)
}
prior_tau_after <- function(x){
  dexp(x,rate=10.0,log=TRUE)
}
priorsEpisodicBD <- c("diversification_before"=prior_delta_before, "turnover_before"=prior_tau_before, "diversification_after"=prior_delta_after, "turnover_after"=prior_tau_after)

rateChangeTime <- max(times)/2

episodicBD <- function(params){
  speciation <- c(params[1]+params[2], params[3]+params[4])
  extinction <- c(params[2], params[4])
  lnl <- tess.likelihood.rateshift(times, lambda=speciation, mu=extinction, rateChangeTimesLambda=rateChangeTime, rateChangeTimesMu=rateChangeTime, samplingProbability=sampfrac, log=TRUE)
  return(lnl)
}

mlepisodicBD <- tess.steppingStoneSampling(likelihoodFunction = episodicBD, priors = priorsEpisodicBD, parameters = runif(4,0,1), logTransforms = c(TRUE,TRUE,TRUE,TRUE), iterations = 1000, burnin = 100, K = 50)

#Mass extinction birth-death on S&B tree (see Hohna TESS vignette)
prior_delta <- function(x){
  dexp(x, rate=10.0, log=TRUE)
}
prior_tau <- function(x){
  dexp(x, rate=10.0, log=TRUE)
}
prior_time <- function(x){
  dunif(x, min=max(times)/2, max=max(times), log=TRUE)
}
priorsMassExtinctionBD <- c("diversification"=prior_delta, "turnover"=prior_tau, "mass-extinctiontime"=prior_time)

survivalProbability <- 0.1

likelihoodMassExtinctionBD <- function(params) { 
  speciation <- params[1] + params[2] 
  extinction <- params[2] 
  time <- params[3] 
  lnl <- tess.likelihood(times, lambda=speciation, mu=extinction, massExtinctionTimes=time, massExtinctionSurvivalProbabilities=survivalProbability, samplingProbability=1.0, log=TRUE) 
  return(lnl) 
}

mlMassExtinctionBD <- tess.steppingStoneSampling(likelihoodFunction=likelihoodMassExtinctionBD, priors=priorsMassExtinctionBD, parameters=c(runif(2,0,1), max(times)*3/4), logTransforms=c(TRUE,TRUE,FALSE), iterations=1000, burnin=100, K=50)

#Compare candidate models
candidateModels <- c("ConstBD"=mlconstBD, "EpisodicBD"=mlepisodicBD, "MassExtBD"=mlMassExtinctionBD)
marginalLikelihoodGrid <- expand.grid(M0=names(candidateModels), M1=names(candidateModels))
marginalLikelihoodGrid$BF <- 2 * (candidateModels[marginalLikelihoodGrid$M0] - candidateModels[marginalLikelihoodGrid$M1])
marginalLikelihoodGrid <- marginalLikelihoodGrid[order(marginalLikelihoodGrid$BF, decreasing=TRUE),]
print(marginalLikelihoodGrid)

#########################################################################################################
#Can we reproduce their results with a tree that we have removed taxa from?
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


