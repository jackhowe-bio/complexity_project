# Model 6: germline vs number of cell types / phylogeny WITHOUT CELL NUMBER
#This is the same as model4, but with the number of cells as explanatory variable removed
#to test whether the cell number is important

setwd('/Users/pcx971/Documents/oxford/complexity/complexity_project/HPC_Analyses_Genus/')

## load packages
library(ape)
library(tidyverse)
library(MCMCglmm)
library(parallel)
library(coda)
library(ggthemes)


#import the MCMCglmm parameters
iterations<- readRDS('RScripts/R_Objects/iterations.RDS')
burnin<- readRDS('RScripts/R_Objects/burnin.RDS')
thinning <- readRDS('RScripts/R_Objects/thinning.RDS')
n_chains <- readRDS('RScripts/R_Objects/n_chains.RDS')
################################

#Read in the data and the tree
df = readRDS('RScripts/R_Objects/metadata.RDS')
inv_tree = readRDS('RScripts/R_Objects/inv_tree.RDS')
df<- subset(df, !is.na(GermTimeSimp ))

#Read in the priors
p1=readRDS('RScripts/R_Objects/priors_1.RDS') #sets prior for residual variance, the defaults are used as priors for fixed effects (see MCMCglmm course notes)
p2=readRDS('RScripts/R_Objects/priors_2.RDS') #sets prior for residual variance, the defaults are used as priors for fixed effects (see MCMCglmm course notes)
p3=readRDS('RScripts/R_Objects/priors_3.RDS') #sets prior for residual variance, the defaults are used as priors for fixed effects (see MCMCglmm course notes)


i = 1
for(prior_set in list(p1,p2,p3)){
  Model4_ROTL_parallel <- mclapply(1:n_chains, function(i){
    ## EDITED HERE: Removed cell number as explanatory variable
    MCMCglmm(Types ~ GermTimeSimp-1, #-1 here removes the intercept equivalent to 0 in brms
             random = ~species_rotl, ginverse=list(species_rotl=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
             family ="poisson",data = df,prior=prior_set, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=T)
  }, mc.cores = 6)
  names(Model4_ROTL_parallel)<- c('chain1','chain2','chain3','chain4', 'chain5','chain6')
  
  
  Model4_ROTL_Sol<- mcmc.list(lapply(Model4_ROTL_parallel, function(m) m$Sol))
  #plot(Model4_ROTL_Sol) 
  gelman.diag(Model4_ROTL_Sol,multivariate = FALSE)
  
  Model4_ROTL_VCV<- mcmc.list(lapply(Model4_ROTL_parallel, function(m) m$VCV))
  #plot(Model4_ROTL_VCV) 
  gelman.diag(Model4_ROTL_VCV,multivariate = FALSE)
  chain1<- Model4_ROTL_parallel$chain1
  summary(Model4_ROTL_parallel$chain1)
  
  M4_Sol_gg<- ggmcmc::ggs(Model4_ROTL_Sol)
  
  diagnostic_filename<- paste('RScripts/ModelDiagnostics/p', i,'/Model4_ROTL_p',i, '.pdf', sep = '')
  #ggmcmc::ggmcmc(M4_Sol_gg, file = diagnostic_filename)
  
  MCMCglmm_filename<- paste('RScripts/ModelOutputs/p', i,'/Model4_ROTL_p',i, '.RDS', sep = '')
  saveRDS(Model4_ROTL_parallel, MCMCglmm_filename)
  i = i + 1
}