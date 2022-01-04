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
M2.1 <- readRDS('RScripts/model_outputs/ModelMulti.RDS')

#import the MCMCglmm parameters
iterations<- 20000000
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
M2.1<- mclapply(1:n_chains, function(i){
  MCMCglmm(cbind(FissionBinary,EarlyGermlineBinary) ~ trait-1,
                random = ~us(trait):species, #2x2 phylogenetic covariance matrix
                rcov = ~us(trait):units, #2x2 residual covariance matrix
                ginverse=list(species=inv_tree),family = c("categorical","categorical"), data = df_binary,prior=pM2.1, nitt=12000000, burnin=burnin, thin=thinning,verbose = T)
}, mc.cores = 4)

#Univariate regression is covariance between response and predictor/variance in predictor
M2.1_VCV<- mcmc.list(lapply(M2.1, function(m) m$VCV))
plot(M2.1_VCV) 
gelman.diag(M2.1_VCV,multivariate = FALSE)

M2.1_VCV_gg<- ggmcmc::ggs(M2.1_VCV)

# for between species use the parameters ending in species 
Covariance<- M2.1_VCV_gg %>%
  filter(Parameter == 'traitFissionBinary.1:traitEarlyGermlineBinary.1.species') %>% 
  rename(Covariance = value) %>% select(-Parameter)

VariancePredictor<- M2.1_VCV_gg %>%
  filter(Parameter == 'traitFissionBinary.1:traitFissionBinary.1.species') %>%
  rename(PredictorVariance = value) %>% select(-Parameter)

Regression<- merge(Covariance, VariancePredictor) %>%
  mutate(value = Covariance / PredictorVariance)
glimpse(Regression)
ggplot(Regression, aes(x = Iteration, y = value, colour = as.factor(Chain))) + geom_line(size = 1, alpha = 0.5) + geom_smooth() + theme_minimal()

# units?

saveRDS(M2.1, 'RScripts/model_outputs/ModelMulti.RDS')

