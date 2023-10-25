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

############################################################################################
#Extract ages of all the nodes in the phylogeny
node_ages <- dateNodes(sb_tree_trimmed)
node_ages <- sort(node_ages, decreasing = T)
node_ages[1] #This is the root, dated to approximately 135.9 million years ago
length(which(node_ages > 65)) / (length(node_ages) / 2) #Only 0.35% of all evolutionary events occur at or before the (generous) 65mya mark

#How many lineages exist at 65 mya?
cretaceous_slice <- timeSliceTree(ttree=sb_tree_trimmed, sliceTime=65, plot=T)
length(cretaceous_slice$tip.label) / length(sb_tree_trimmed$tip.label)

