# Ancestral State Reconstructions using Charlie's Bayes functions
setwd('/Users/pcx971/Documents/oxford/complexity/complexity_project/HPC_Analyses/')


# load the files containing the functions
devtools::source_url("https://raw.githubusercontent.com/charliecornwallis/Rfunctions/master/ASR_functions.R")
source('RScripts/pMCMC_Script.R')

#load packages
library(ape)
library(MCMCglmm)
library(tidyr)
library(parallel)

# Read in the data and the tree
df = readRDS('RScripts/R_Objects/metadata.RDS')
inv_tree = readRDS('RScripts/R_Objects/inv_tree.RDS')
tree<- read.tree('anvio/Phylogeny_ROTL.txt')

#import the MCMCglmm parameters
iterations<- readRDS('RScripts/R_Objects/iterations.RDS')
burnin<- readRDS('RScripts/R_Objects/burnin.RDS')
thinning <- readRDS('RScripts/R_Objects/thinning.RDS')
n_chains <- readRDS('RScripts/R_Objects/n_chains.RDS')
################################


#import the MCMCglmm parameters
iterations<- 10000
burnin<- 1000
thinning <- 10
n_chains <- 3
################################

# Running models to reconstruct the traits

# cell types
#Read in the priors
p1=readRDS('RScripts/R_Objects/priors_1.RDS') #sets prior for residual variance, the defaults are used as priors for fixed effects (see MCMCglmm course notes)
p2=readRDS('RScripts/R_Objects/priors_2.RDS') #sets prior for residual variance, the defaults are used as priors for fixed effects (see MCMCglmm course notes)
p3=readRDS('RScripts/R_Objects/priors_3.RDS') #sets prior for residual variance, the defaults are used as priors for fixed effects (see MCMCglmm course notes)

priors= list(p1,p2,p3)
priors= list(p1)


i = 1
for(prior in priors){
  # run MCMCglmm
  print(prior)
  MCMCglmm_Types <- mclapply(1:n_chains, function(i){
    MCMCglmm(Types ~ 1, #-1 here removes the intercept equivalent to 0 in brms
             random = ~species_rotl, ginverse=list(species_rotl=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
             family ="poisson",data = df,prior=prior, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=T, pl = T)
  }, mc.cores = 3)
  
  names(MCMCglmm_Types)<- c('chain1','chain2','chain3')

  MCMCglmm_filename<- paste('AncestralStates_p', i,'_CellTypes.RDS', sep = '')
  saveRDS(MCMCglmm_Types, MCMCglmm_filename)
  i = i + 1
  


MCMCglmm_Types<- readRDS(MCMCglmm_filename)


predTypes<-pred_states_mcmcglmm(mr="No", # is the model multi-response?
                           trees=tree, # what is the phylogeny used?
                           phy_name="species_rotl", # the in-model name for the phylogeny?
                           model=MCMCglmm_Types$chain2, # the model name
                           dat=as.data.frame(df), #the dataframe
                           trait1="Type", #the response variable
                           binomial="No") # whether or not this trait is binary -> since I’m reconstructing our continuous response variable, set to ‘no’
#state1="early", # name of first state of binary trait
#state2="late", # name of second state of trait
#cutoff=0.5)

print(predTypes)

}

# cell numbers
i = 1
for(prior_set in priors){
  # run MCMCglmm
  MCMCglmm_Number <- mclapply(1:n_chains, function(i){
    MCMCglmm(Number ~ 1, #-1 here removes the intercept equivalent to 0 in brms
             random = ~species_rotl, ginverse=list(species_rotl=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
             family ="poisson",data = df,prior=prior_set, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=T, pl = T)
  }, mc.cores = 3)
  
  names(MCMCglmm_Number)<- c('chain1','chain2','chain3')
  
  MCMCglmm_filename<- paste('AncestralStates_p', i,'_CellNumbers.RDS', sep = '')
  saveRDS(MCMCglmm_Number, MCMCglmm_filename)
  i = i + 1
  


MCMCglmm_Number<- readRDS(MCMCglmm_filename)


predNumber<-pred_states_mcmcglmm(mr="No", # is the model multi-response?
                                trees=tree, # what is the phylogeny used?
                                phy_name="species_rotl", # the in-model name for the phylogeny?
                                model=MCMCglmm_Number$chain1, # the model name
                                dat=as.data.frame(df), #the dataframe
                                trait1="Number", #the response variable
                                binomial="No") # whether or not this trait is binary -> since I’m reconstructing our continuous response variable, set to ‘no’
#state1="early", # name of first state of binary trait
#state2="late", # name of second state of trait
#cutoff=0.5)
predNumber$`log(estimate)`<- log(predNumber$estimate)
print(predNumber)
}


# fission
i = 1
for(prior_set in priors){
  # run MCMCglmm
  MCMCglmm_Fission <- mclapply(1:n_chains, function(i){
    MCMCglmm(as.factor(Fission) ~ 1, #-1 here removes the intercept equivalent to 0 in brms
             random = ~species_rotl, ginverse=list(species_rotl=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
             family ="categorical",data = df,prior=prior_set, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=T, pl = T)
  }, mc.cores = 3)
  
  names(MCMCglmm_Fission)<- c('chain1','chain2','chain3')
  
  MCMCglmm_filename<- paste('AncestralStates_p', i,'_Fission.RDS', sep = '')
  saveRDS(MCMCglmm_Fission, MCMCglmm_filename)
  i = i + 1
  


MCMCglmm_Fission<- readRDS(MCMCglmm_filename)


predFission<-pred_states_mcmcglmm(mr="No", # is the model multi-response?
                                 trees=tree, # what is the phylogeny used?
                                 phy_name="species_rotl", # the in-model name for the phylogeny?
                                 model=MCMCglmm_Fission$chain1, # the model name
                                 dat=as.data.frame(df), #the dataframe
                                 trait1="Fission", #the response variable
                                 binomial="Yes", # whether or not this trait is binary -> since I’m reconstructing our continuous response variable, set to ‘no’
                                 state1="0", # name of first state of binary trait
                                 state2="1", # name of second state of trait
                                cutoff=0.2)

print(predFission)
}

# germlines
i = 1
for(prior_set in priors){
  # run MCMCglmm
  MCMCglmm_Germline <- mclapply(1:n_chains, function(i){
    MCMCglmm(GermNumeric ~ 1, #-1 here removes the intercept equivalent to 0 in brms
             random = ~species_rotl, ginverse=list(species_rotl=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
             family ="categorical",data = df,prior=p1, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=T, pl = T)
  }, mc.cores = 3)
  
  names(MCMCglmm_Germline)<- c('chain1','chain2','chain3')
  
  MCMCglmm_filename<- paste('AncestralStates_p', i,'_Germline.RDS', sep = '')
  saveRDS(MCMCglmm_Germline, MCMCglmm_filename)
  i = i + 1
  
  
  
  MCMCglmm_Germline<- readRDS(MCMCglmm_filename)
  
  
  predFission<-pred_states_mcmcglmm(mr="No", # is the model multi-response?
                                    trees=tree, # what is the phylogeny used?
                                    phy_name="species_rotl", # the in-model name for the phylogeny?
                                    model=MCMCglmm_Germline$chain1, # the model name
                                    dat=as.data.frame(df), #the dataframe
                                    trait1="GermNumeric", #the response variable
                                    binomial="yes", # whether or not this trait is binary -> since I’m reconstructing our continuous response variable, set to ‘no’
                                    state1="1", # name of first state of binary trait
                                    state2="2", # name of second state of trait
                                    cutoff=0.5)
  
  predFission
}

