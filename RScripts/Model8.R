# Model 8: germline vs number of cells w/ phylogeny

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


#Setting the priors (there are others to try especially for binary data, but these usually work well for other familys)
p2=list(R = list(V = 1, nu=0.002), 
        G = list(G1=list(V=1, nu=0.002)))#prior for random variances - there is G1 because only 1 random effect. If there were 2 it would be G = list(G1=list(V=1, nu=0.002),G2=list(V=1, nu=0.002))


Model8_parallel <- mclapply(1:n_chains, function(i){
  MCMCglmm(cell_number ~ cell_types ~ germline_timing_simple-1 + scale(log(cell_number)), #-1 here removes the intercept equivalent to 0 in brms
           random = ~species, ginverse=list(species=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
           family ="poisson",data = df,prior=p2, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=T)
}, mc.cores = 4)
names(Model8_parallel)<- c('chain1','chain2','chain3','chain4', 'chain 5','chain 8')


Model8_Sol<- mcmc.list(lapply(Model8_parallel, function(m) m$Sol))
#plot(Model8_Sol) 
gelman.diag(Model8_Sol,multivariate = FALSE)

Model8_VCV<- mcmc.list(lapply(Model8_parallel, function(m) m$VCV))
#plot(Model8_VCV) 
gelman.diag(Model8_VCV,multivariate = FALSE)
chain1<- Model8_parallel$chain1
summary(Model8_parallel$chain1)

M8_Sol_gg<- ggmcmc::ggs(Model8_Sol)
ggmcmc::ggmcmc(M8_Sol_gg, file = 'MCMCglmmDiagnostics/Model8.pdf')

summary(Model8_parallel$chain1)
saveRDS(Model8_parallel, 'RScripts/model_outputs/Model8.RDS')
####################

####SAVE OUTPUT #### 
#save image of environment so it doesn't need run each time
#save.image(file = 'RScripts/Model8.Rdata')