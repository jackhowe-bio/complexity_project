# Model 5: fission vs number of cells / phylogeny

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


Model5_parallel <- mclapply(1:n_chains, function(i){
  MCMCglmm(cell_number ~FissionOrBuddingObserved_Species-1, #-1 here removes the intercept equivalent to 0 in brms
           random = ~species, ginverse=list(species=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
           family ="poisson",data = df,prior=p2, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=T)
}, mc.cores = 4)
names(Model5_parallel)<- c('chain1','chain2','chain3','chain4', 'chain 5','chain 6')


Model5_Sol<- mcmc.list(lapply(Model5_parallel, function(m) m$Sol))
#plot(Model5_Sol) 
gelman.diag(Model5_Sol,multivariate = FALSE)

Model5_VCV<- mcmc.list(lapply(Model5_parallel, function(m) m$VCV))
#plot(Model5_VCV) 
gelman.diag(Model5_VCV,multivariate = FALSE)
chain1<- Model5_parallel$chain1
summary(Model5_parallel$chain1)

M5_Sol_gg<- ggmcmc::ggs(Model5_Sol)
ggmcmc::ggmcmc(M5_Sol_gg, file = 'MCMCglmmDiagnostics/Model5.pdf')

summary(Model5_parallel$chain1)
saveRDS(Model5_parallel, 'RScripts/model_outputs/Model5.RDS')
####################

####SAVE OUTPUT #### 
#save image of environment so it doesn't need run each time
#save.image(file = 'RScripts/Model5.Rdata')