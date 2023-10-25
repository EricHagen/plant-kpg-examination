#######################################################################################################
# How many fossils does it take to get an accurate TESS result?
#######################################################################################################

#Load necessary packages
x <- c("diversitree", "geiger", "stringr", "castor")
invisible(lapply(x, library, character.only=T))

#Set your directory & source some functions:
folder <- "/Users/erichagen1/Desktop/University_of_Toronto_Stuff/publications/plant-kpg-examination/"
source(paste0(folder, "kpg_simulations.R"))

#IF YOU ARE RUNNING THIS THE FIRST TIME
extinct_removed_orig <- NULL
while(is.null(extinct_removed_orig) || length(extinct_removed_orig$tip.label) < 100){
  extinct_removed_orig <- tree.bd(pars=c(1, 0.3), max.taxa=10000, include.extinct=TRUE)
}
save(extinct_removed_orig, file=paste0(folder, "extinct_removed_orig.RData"))

#OTHERWISE, SAVE THIS TREE AND LOAD IT:
load(paste0(folder, "extinct_removed_orig.RData"))

#Dropping ALL extinct taxa from the tree:
extinct_removed <- drop.extinct(extinct_removed_orig)

#Removing a certain number of extinct taxa:
extinct_taxa <- is.extinct(extinct_removed_orig)
howmany = 1
extinct_taxa_toremove <- extinct_taxa[1:(length(extinct_taxa - howmany))]
mixed_phy <- drop.tip(extinct_removed_orig, which(!extinct_removed_orig$tip.label %in% extinct_taxa_toremove))

#Or, remove taxa at random:
extinct_taxa_toremove <- sample(extinct_taxa, howmany)
mixed_phy <- drop.tip(extinct_removed_orig, extinct_taxa_toremove)

#Or, remove taxa by clade:
extinct_taxa_toremove <- vector(mode="character", length=howmany)
removed <- sample(extinct_taxa, 1)
possibles <- nodepath(extinct_removed_orig, find_root(extinct_removed_orig), as.numeric(gsub("ex", "", removed)))
clade_sizes <- vector(mode="numeric", length=length(possibles))
for(i in 1:length(possibles)){
  descs <- geiger::tips(extinct_removed_orig, possibles[i])
  extincts <- str_detect(descs, "ex")
  clade_sizes[i] <- length(descs[which(extincts==TRUE)])
}
clade_sizes <- clade_sizes[which(clade_sizes > howmany)]
removeable_tips <- geiger::tips(extinct_removed_orig, possibles[length(clade_sizes)])
extincts <- str_detect(removeable_tips, "ex")
extinct_taxa_toremove <- sample(removeable_tips[which(extincts==TRUE)], howmany)
mixed_phy <- drop.tip(extinct_removed_orig, extinct_taxa_toremove)

#Run the TESS model fitting function
grid_extant <- model_comparison(sampfrac = (length(extinct_removed$tip.label) / length(extinct_removed_orig$tip.label)), phy = extinct_removed)
grid_mixed <- model_comparison(sampfrac = (length(mixed_phy$tip.label) / length(extinct_removed_orig$tip.label)), phy = mixed_phy)

#How do Bayes factors compare in both scenarios?
colnames(grid_extant) <- c("M0_All_Extant", "M1_All_Extant", "BF_All_Extant")
colnames(grid_mixed) <- c("M0_Mixed", "M1_Mixed", "BF_Mixed")
grid_all <- cbind(grid_extant, grid_mixed)
grid_all

#Do this in a loop


