# Script that fits MCMCglmm models with varying combinations or parameters for 

# optimising analyses
setwd('/Users/pcx971/Documents/oxford/complexity/complexity_project/HPC_Analyses_Genus/')
library(ape)
library(tidyverse)
library(MCMCglmm)
library(parallel)
library(coda)
library(ggthemes)
library(ggpubr)
        
#Read in the data and the tree
df = readRDS('RScripts/R_Objects/metadata.RDS')
inv_tree = readRDS('RScripts/R_Objects/inv_tree.RDS')

#Read in the priors
p1=readRDS('RScripts/R_Objects/priors_1.RDS') #sets prior for residual variance, the defaults are used as priors for fixed effects (see MCMCglmm course notes)

# define the testing arguments by generating dataframe
testing_args <- expand.grid(n_itts = as.integer(seq(500000, 10000000, length = 7)), n_thin = as.integer(seq(100,1000, length = 2)), n_burn = as.integer(seq(100000, 1000000, length = 3)))
#make sure the burn in length does not exceed the total length (otherwise error)
testing_args<- subset(testing_args, n_itts >= 2*n_burn)

#running MCMCglmm for all combinations
OptimisationRunParallel <- mclapply(1:nrow(testing_args), function(i, n_itts,n_burn, n_thin){
  MCMCglmm(Number ~Fission-1, #-1 here removes the intercept equivalent to 0 in brms
           random = ~species_rotl, ginverse=list(species_rotl=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
           family ="poisson",data = df ,prior=p1, nitt=n_itts[i], burnin=n_burn[i], thin=n_thin[i], verbose = F, pr=T)
}, mc.set.seed=T, mc.cores = 6, n_itts = testing_args$n_itts, n_burn = testing_args$n_burn, n_thin = testing_args$n_thin)

options(scipen = 8) #makes sure that numbers come out properly, penalty for scientific notiation under 8 digits
names(OptimisationRunParallel)<- paste('it:',testing_args$n_itts,'thin:',testing_args$n_thin,'burnin:',testing_args$n_burn, sep = '') #set a name for each run that gives the parameter info for later
summary(OptimisationRunParallel) #summary of run

autocorrelation_tests_VCV<- data.frame(sapply(OptimisationRunParallel, function(m) autocorr(m$VCV, lags = 1))) #conduct autocorrelation tests for variance

#make the data readable for ggplot
autocorrelation_tests_VCV_long<- autocorrelation_tests_VCV %>% 
  mutate(info = rownames(.)) %>%
  extract(info, 
          into = c('iterations','thinning','burnin'),
          regex = "it:(.*)thin:(.*)burnin:(.*)") %>%
  rename(autocorrelation = 1) %>%
  mutate(iterations = as.numeric(iterations))


autocorrelation_tests_Sol<- data.frame(sapply(OptimisationRunParallel, function(m) autocorr(m$Sol, lags = 1))) #conduct autocorrelation tests for mean
autocorrelation_tests_Sol_long<- autocorrelation_tests_Sol %>% 
  pivot_longer(everything(), 
               names_pattern = "it\\.(.*)thin\\.(.*)burnin\\.(.*)",
               names_to = c('iterations','thinning','burnin'),
               values_to = "autocorrelation") %>%
  mutate(iterations = as.numeric(iterations))

#plot the figures
OptimisingVariancePlot<- ggplot(autocorrelation_tests_VCV_long, aes(x = iterations, y = abs(autocorrelation), colour = burnin, shape = thinning)) + 
  geom_point() + theme_minimal() + labs(title = "Variance", x = "Iterations", y = "Absolute autocorrelation" )
  
OptimisingMeanPlot<- ggplot(autocorrelation_tests_Sol_long, aes(x = iterations, y = abs(autocorrelation), colour = burnin, shape = thinning)) + 
  geom_point() + theme_minimal() + labs(title = "Mean", x = "Iterations", y = "Absolute autocorrelation" )


OptimisingPlot<- ggarrange(OptimisingVariancePlot,OptimisingMeanPlot, ncol = 1, nrow = 2)

pdf('figures/Optimisation.pdf')
OptimisingPlot
dev.off()

#save image of environment so it doesn't need run each time
saveRDS(OptimisationRunParallel, 'RScripts/R_Objects/OptimisingModels.RDS')