# Model 6: Fission vs number of cell types w/ phylogeny

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


Model6_parallel <- mclapply(1:n_chains, function(i){
  MCMCglmm(cell_types ~ FissionOrBuddingObserved_Species-1 + scale(log(cell_number)), #-1 here removes the intercept equivalent to 0 in brms
           random = ~species, ginverse=list(species=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
           family ="poisson",data = df,prior=p2, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=T)
}, mc.cores = 4)
names(Model6_parallel)<- c('chain1','chain2','chain3','chain4', 'chain 5','chain 6')


Model6_Sol<- mcmc.list(lapply(Model6_parallel, function(m) m$Sol))
#plot(Model6_Sol) 
gelman.diag(Model6_Sol,multivariate = FALSE)

Model6_VCV<- mcmc.list(lapply(Model6_parallel, function(m) m$VCV))
#plot(Model6_VCV) 
gelman.diag(Model6_VCV,multivariate = FALSE)
chain1<- Model6_parallel$chain1
summary(Model6_parallel$chain1)

M6_Sol_gg<- ggmcmc::ggs(Model6_Sol)
ggmcmc::ggmcmc(M6_Sol_gg, file = 'MCMCglmmDiagnostics/Model6.pdf')

summary(Model6_parallel$chain1)
saveRDS(Model6_parallel, 'RScripts/model_outputs/Model6.RDS')
####################

####SAVE OUTPUT #### 
#save image of environment so it doesn't need run each time
#save.image(file = 'RScripts/Model6.Rdata')