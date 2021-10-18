# Model 3: Germline vs number of cells

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

M3_parallel <- mclapply(1:n_chains, function(i){
  MCMCglmm(cell_number ~ germline_timing_simple-1, #-1 here removes the intercept equivalent to 0 in brms
           family ="poisson",data = df,prior=p1, nitt=iterations, burnin=burnin, thin=thinning,verbose = F,pr=T)
}, mc.cores = 4)
names(M3_parallel)<- c('chain1','chain2','chain3','chain4','chain5','chain6')

M3_Sol<- mcmc.list(lapply(M3_parallel, function(m) m$Sol))
plot(M3_Sol) 
gelman.diag(M3_Sol,multivariate = FALSE)

M3_VCV<- mcmc.list(lapply(M3_parallel, function(m) m$VCV))
plot(M3_VCV) 
gelman.diag(M3_VCV,multivariate = FALSE)

M3_Sol_gg<- ggmcmc::ggs(M3_Sol)
ggmcmc::ggmcmc(M3_Sol_gg, file = 'MCMCglmmDiagnostics/Model3.pdf')



summary(M3_parallel$chain1)
saveRDS(M3_parallel, 'RScripts/model_outputs/Model3.RDS')
####################

####SAVE OUTPUT #### 
#save image of environment so it doesn't need run each time
#save.image(file = 'RScripts/Model3.Rdata')