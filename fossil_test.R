#######################################################################################################
# How many fossils does it take to get an accurate TESS result?
#######################################################################################################

#Load necessary packages
x <- c("diversitree", "geiger", "stringr", "castor", "parallel", "dplyr", "rlist")
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

#Ladderize the tree:
extinct_removed <- ladderize(extinct_removed)

#Removing a certain number of extinct taxa:
extinct_taxa <- is.extinct(extinct_removed_orig)
howmany = 1
extinct_taxa_toremove <- extinct_taxa[1:(length(extinct_taxa - howmany))]
mixed_phy <- drop.tip(extinct_removed_orig, which(!extinct_removed_orig$tip.label %in% extinct_taxa_toremove))
mixed_phy <- ladderize(mixed_phy)

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

#Other trees:
mass_ext_phy <- NULL
while(is.null(mass_ext_phy)){
  mass_ext_phy <- sim.rateshift.taxa(n=5000, numbsim=1, lambda=c(1,0.5), mu=c(0.3, 0.4), frac=c(1,0.1), times=c(0,30))
}
mass_ext_extant <- drop.extinct(mass_ext_phy[[1]])

#######################################################################################################
# Do this in a loop
#######################################################################################################

howmany = 13
numz <- seq(10001, length(extinct_removed_orig$tip.label), by=howmany)
numz <- numz[-length(numz)]
iterationz <- seq(1, length(numz), by=1)
all_results <- vector(mode="list", length=length(numz))
extinct_taxa <- is.extinct(extinct_removed_orig)
extinct_taxa_toremove <- sample(extinct_taxa, (length(extinct_removed_orig$tip.label)-10001))
#for(i in 1:length(numz)){
#while(length(mixed_phy$tip.label) > (extinct_removed$tip.label + 1)){
get_grids <- function(i, numz, howmany, extinct_taxa_toremove, extinct_removed_orig){
  j=numz[i]
  mixed_phy <- extinct_removed_orig
  tipz <- extinct_taxa_toremove[1:(j*howmany)]
  mixed_phy <- drop.tip(mixed_phy, tipz)
  grid_mixed <- model_comparison(sampfrac = (length(mixed_phy$tip.label) / length(extinct_removed_orig$tip.label)), phy = mixed_phy)
  grid_mixed$Taxa_Removed <- rep((length(extinct_removed_orig$tip.label) - length(mixed_phy$tip.label)), nrow(grid_mixed))
  #all_results <- rbind(all_results, grid_mixed)
  write.csv(grid_mixed, file=paste0("/home/erik/Documents/tess_kpg_project/tess_runs/grid", i, ".csv"))
}

#Parallelize the loop
ncores = 13
mclapply(iterationz, function(x) get_grids(i=x, numz=numz, howmany=howmany, extinct_taxa_toremove=extinct_taxa_toremove, extinct_removed_orig=extinct_removed_orig), mc.cores=ncores)
#mclapply(iterationz, function(x) get_grids(i=x, numz=numz, howmany=howmany, extinct_taxa_toremove=extinct_taxa_toremove, extinct_removed_orig=mass_ext_phy[[1]]), mc.cores=ncores)

#Convert list to dataframe
all_results_df <- as.data.frame(matrix(nrow=0, ncol=4))
colnames(all_results_df) <- c("M0", "M1", "BF", "Taxa_Removed")
for(i in iterationz){
  curr <- read.csv(paste0("/home/erik/Documents/tess_kpg_project/tess_runs/grid", i, ".csv"))
  curr <- curr[,c(2:5)]
  all_results_df <- rbind(all_results_df, curr)
}

#Examine results
results_df_bd <- read.csv("/home/erik/Documents/tess_kpg_project/tess_runs/results_files/birth_death_fossiladdtest.csv")
results_df_me <- read.csv("/home/erik/Documents/tess_kpg_project/tess_runs/results_files/mass_ext_fossiladdtest.csv")
results_df_bd <- results_df_bd[,c(2:5)]
results_df_me <- results_df_me[,c(2:5)]

#Altering each dataframe because some stuff messed up

#Keeping in mind that:
#(1) Each TESS run produces 9 rows in each dataframe
#(2) The mass_ext tree has 7257 tips total and 5000 tips without fossils (2256 between)
#(3) The bd tree has 14330 tips total and 10000 tips without fossils (4329 between)
#(4) The mass_ext tree removed 16 tips at a time, while the bd tree removed 13 tips at a time. Let's fix
#the "Taxa_Removed" columns in each results dataframe:
fix_seq <- seq(14330, 10001, by=-13)
fix_seq <- fix_seq - 10001
fix_seq <- fix_seq[-length(fix_seq)]
results_df_bd$Taxa_Removed <- rep(fix_seq, each=9)
fix_seq <- seq(7257, 5001, by=-16)
fix_seq <- fix_seq - 5001
fix_seq <- fix_seq[-length(fix_seq)]
results_df_me$Taxa_Removed <- rep(fix_seq, each=9)
#Now we make separate dataframes to check the accuracy of each run
bd_check <- as.data.frame(matrix(ncol=3, nrow=333))
me_check <- as.data.frame(matrix(ncol=3, nrow=141))
colnames(bd_check) <- colnames(me_check) <- c("Favored_Model", "BF", "Taxa_Removed")
bd_seq <- seq(from=1, to=2989, by=9)
for(i in 1:length(bd_seq)){
  j = bd_seq[i]
  bd_check[i,] <- results_df_bd[j, c(2:4)]
}
#Constant birth-death consistently favored across all birth-death trees, regardless of taxa removed
me_seq <- seq(from=1, to=1261, by=9)
for(i in 1:length(me_seq)){
  j = me_seq[i]
  me_check[i,] <- results_df_me[j, c(2:4)]
}
##Constant birth-death consistently favored across all mass_ext trees(!!!!), regardless of taxa removed

#Bayes Factors do not change very much at all for either tree as taxa are added back

