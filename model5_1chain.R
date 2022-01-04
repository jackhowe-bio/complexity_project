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


M5<- MCMCglmm(cell_number ~FissionOrBuddingObserved_Species-1, #-1 here removes the intercept equivalent to 0 in brms
           random = ~species, ginverse=list(species=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
           family ="poisson",data = df,prior=p2, nitt=500000, burnin=50000, thin=10 ,verbose = T, pr=T)

pacman::p_load(devtools)
#source_url("https://raw.githubusercontent.com/charliecornwallis/Rfunctions/master/MCMCglmmProc.R")
summary(M5)
source('RScripts/pMCMCglmmScript.R')

fixed_del_phy<- colnames(M5$Sol)[grep('species.', colnames(M5$Sol))]

#Apply function (description of function terms can be found on github: https://github.com/charliecornwallis/Rfunctions/blob/master/MCMCglmmProc.R)
SItablesXL<-MCMCglmmProc(model=M5,link="poisson",start_row=1,create_sheet="yes",sheet="Table S1",title="Table S1: Jack project", 
                         fixed_names=c("FissionOrBuddingObserved_Species0", "FissionOrBuddingObserved_Species1","FissionOrBuddingObserved_Species?", "FissionOrBuddingObserved_Species"),
                         fixed_diffinc = c("all"),
                         Include_random = "yes",
                         pvalues = "include",
#                         randomvar_names=c("species","units"),
                         fixed_del=c(fixed_del_phy),
                         variances=NULL,padding=3)




####################

####SAVE OUTPUT #### 
#save image of environment so it doesn't need run each time
#save.image(file = 'RScripts/Model5.Rdata')