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
M2.1<- mclapply(1:6, function(i){
  MCMCglmm(cbind(FissionBinary,EarlyGermlineBinary) ~ trait-1,
                random = ~us(trait):species, #2x2 phylogenetic covariance matrix
                rcov = ~us(trait):units, #2x2 residual covariance matrix
                ginverse=list(species=inv_tree),family = c("categorical","categorical"), data = df_binary,prior=pM2.1, nitt=600000, burnin=10000, thin=100,verbose = T)
}, mc.cores = 4)

M2.1_Sol<- mcmc.list(lapply(M2.1, function(m) m$Sol))
plot(M2.1_Sol) 
gelman.diag(M2.1_Sol,multivariate = FALSE)

M2.1_VCV<- mcmc.list(lapply(M2.1_parallel, function(m) m$VCV))
plot(M2.1_VCV) 
gelman.diag(M2.1_VCV,multivariate = FALSE)

M1_Sol_gg<- ggmcmc::ggs(M2.1_Sol)
ggmcmc::ggmcmc(M2.1_Sol_gg, file = 'MCMCglmmDiagnostics/Model1.pdf')

