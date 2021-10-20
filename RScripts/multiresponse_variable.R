# script to test have multiple response variables

# How to test for phylogenetic correlation between them? 
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

# both responses need to be 0s or 1s for the model to work
df_binary<- df %>%
  mutate(FissionBinary= ifelse(df$FissionOrBuddingObserved_Species == 1, 1, 0), EarlyGermlineBinary= ifelse(df$germline_timing_simple == 'early', 1, 0)) 



#Setting the priors (there are others to try especially for binary data, but these usually work well for other familys)
pM2.1 =list(B=list(mu=c(0,0), V=diag(c(1,1+pi^2/3))),
            R = list(V = diag(2),nu=1, fix=1), #Residual variance not identifiable for binary variables
            G = list(G1=list(V = diag(2), nu = 1, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000)))))#parameter expanded priors, usually good for binary data

#Here both responses need to 0 and 1s or yes and nos
M2.1<- MCMCglmm(cbind(FissionBinary,EarlyGermlineBinary) ~ cell_number-1,
                random = ~us(cell_number):species, #2x2 phylogenetic covariance matrix
                rcov = ~us(cell_number):units, #2x2 residual covariance matrix
                ginverse=list(species=inv_tree),family = c("categorical","categorical"), data = df_binary,prior=pM2.1, nitt=600000, burnin=10000, thin=100,verbose = T)

Model8_parallel <- mclapply(1:2, function(i){
  MCMCglmm(cbind(germline_timing_simple, FissionOrBuddingObserved_Species), #-1 here removes the intercept equivalent to 0 in brms
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