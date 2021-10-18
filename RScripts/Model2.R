# Model 2: Fission vs number of cell

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

M2_parallel <- mclapply(1:n_chains, function(i){
  MCMCglmm(cell_types ~ FissionOrBuddingObserved_Species-1 + scale(log(cell_number)), #-1 here removes the intercept equivalent to 0 in brms
           family ="poisson",data = df,prior=p1, nitt=iterations, burnin=burnin, thin=thinning,verbose = F,pr=T)
}, mc.cores = 4)
names(M2_parallel)<- c('chain1','chain2','chain3','chain4','chain5','chain6')

M2_Sol<- mcmc.list(lapply(M2_parallel, function(m) m$Sol))
plot(M2_Sol) 
gelman.diag(M2_Sol,multivariate = FALSE)

M2_VCV<- mcmc.list(lapply(M2_parallel, function(m) m$VCV))
plot(M2_VCV) 
gelman.diag(M2_VCV,multivariate = FALSE)

M2_Sol_gg<- ggmcmc::ggs(M2_Sol)
ggmcmc::ggmcmc(M2_Sol_gg, file = 'MCMCglmmDiagnostics/Model2.pdf')



summary(M2_parallel$chain1)
saveRDS(M2_parallel, 'RScripts/model_outputs/Model2.RDS')
####################

####SAVE OUTPUT #### 
#save image of environment so it doesn't need run each time
save.image(file = 'RScripts/Model2.Rdata')