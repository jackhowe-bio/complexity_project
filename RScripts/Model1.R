# script that runs the analyses for species level data

## load packages
library(ape)
library(tidyverse)
library(MCMCglmm)
library(parallel)
library(coda)
library(ggthemes)

#run the data import script
source('RScripts/DataImportPhylogenyConstruction.R')

#import the MCMCglmm parameters
iterations<- readRDS('RScripts/iterations.RDS')
burnin<- readRDS('RScripts/burnin.RDS')
thinning <- readRDS('RScripts/thinning.RDS')
n_chains <- readRDS('RScripts/n_chains.RDS')
################################
# Model 1: cell number by fission

#Setting the priors
p1=list(R = list(V = 1, nu=0.002)) #sets prior for residual variance, the defaults are used as priors for fixed effects (see MCMCglmm course notes)

M1_parallel <- mclapply(1:6, function(i){
  MCMCglmm(cell_number ~ FissionOrBuddingObserved_Species-1, #-1 here removes the intercept equivalent to 0 in brms
           family ="poisson",data = df,prior=p1, nitt=iterations, burnin=burnin, thin=thinning,verbose = T,pr=T)
}, mc.cores = 4)
names(M1_parallel)<- c('chain1','chain2','chain3','chain4','chain5','chain6')

M1_Sol<- mcmc.list(lapply(M1_parallel, function(m) m$Sol))
plot(M1_Sol) 
gelman.diag(M1_Sol,multivariate = FALSE)

M1_VCV<- mcmc.list(lapply(M1_parallel, function(m) m$VCV))
plot(M1_VCV) 
gelman.diag(M1_VCV,multivariate = FALSE)

M1_Sol_gg<- ggmcmc::ggs(M1_Sol)
ggmcmc::ggmcmc(M1_Sol_gg, file = 'MCMCglmmDiagnostics/Model1.pdf')



summary(M1_parallel$chain1)
saveRDS(M1_parallel, 'RScripts/model_outputs/model_1.RDS')
####################

####SAVE OUTPUT #### 
#save image of environment so it doesn't need run each time
#save.image(file = 'RScripts/SpeciesModels.Rdata')