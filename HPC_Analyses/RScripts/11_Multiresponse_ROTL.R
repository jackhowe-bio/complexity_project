# Model 5: fission vs number of cells / phylogeny
setwd('/Users/pcx971/Documents/oxford/complexity/complexity_project/HPC_Analyses/')

## load packages
library(ape)
library(tidyverse)
library(MCMCglmm)
library(parallel)
library(coda)
library(ggthemes)


#import the MCMCglmm parameters
iterations<- 16000000
burnin<- readRDS('RScripts/R_Objects/burnin.RDS')
thinning <- readRDS('RScripts/R_Objects/thinning.RDS')
n_chains <- readRDS('RScripts/R_Objects/n_chains.RDS')
################################


#Read in the data and the tree
df = readRDS('RScripts/R_Objects/metadata.RDS')
inv_tree = readRDS('RScripts/R_Objects/inv_tree.RDS')

df_binary<- df %>%
  mutate(FissionBinary= as.factor(ifelse(df$Fission == 1, 1, 0)), EarlyGermlineBinary= as.factor(ifelse(df$GermNumeric == 1, 1, 0)))


#Read in the priors
p1=readRDS('RScripts/R_Objects/priors_1.RDS') #sets prior for residual variance, the defaults are used as priors for fixed effects (see MCMCglmm course notes)
p2=readRDS('RScripts/R_Objects/priors_2.RDS') #sets prior for residual variance, the defaults are used as priors for fixed effects (see MCMCglmm course notes)
p3=readRDS('RScripts/R_Objects/priors_3.RDS') #sets prior for residual variance, the defaults are used as priors for fixed effects (see MCMCglmm course notes)

#Setting the priors (there are others to try especially for binary data, but these usually work well for other familys)
p4=list(B=list(mu=c(0,0), V=diag(c(1,1+pi^2/3))),
            R = list(V = diag(2),nu=1, fix=1), #Residual variance not identifiable for binary variables
            G = list(G1=list(V = diag(2), nu = 1, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000)))))#parameter expanded priors, usually good for binary data

#Here both responses need to 0 and 1s or yes and nos

i = 4
for(prior_set in list(p4)){
# run MCMCglmm
  Model_Correlation<- mclapply(1:n_chains, function(i){
    MCMCglmm(cbind(FissionBinary,EarlyGermlineBinary) ~ trait-1,
             random = ~us(trait):species_rotl, #2x2 phylogenetic covariance matrix
             rcov = ~us(trait):units, #2x2 residual covariance matrix
             ginverse=list(species_rotl=inv_tree),family = c("categorical","categorical"), 
             data = df_binary,prior=p4, nitt=iterations, burnin=burnin, thin=thinning,verbose = T)
  }, mc.cores = 6)
names(Model_Correlation)<- c('chain1','chain2','chain3','chain4', 'chain5','chain6')

Model_Correlation_ROTL_Sol<- mcmc.list(lapply(Model_Correlation, function(m) m$Sol))
#plot(Model_Correlation_ROTL_Sol) 
#gelman.diag(Model_Correlation_ROTL_Sol,multivariate = FALSE)

Model_Correlation_ROTL_VCV<- mcmc.list(lapply(Model_Correlation, function(m) m$VCV))
#plot(Model_Correlation_ROTL_VCV) 
gelman.diag(Model_Correlation_ROTL_VCV,multivariate = FALSE)
chain1<- Model_Correlation$chain1
summary(Model_Correlation$chain1)

Model_Correlation_Sol_gg<- ggmcmc::ggs(Model_Correlation_ROTL_Sol)

diagnostic_filename<- paste('RScripts/ModelDiagnostics/p', i,'/Model_Correlation_ROTL_p',i, '.pdf', sep = '')
ggmcmc::ggmcmc(Model_Correlation_Sol_gg, filename = diagnostic_filename)

MCMCglmm_filename<- paste('RScripts/ModelOutputs/p', i,'/Model_Correlation_ROTL_p',i, '.RDS', sep = '')
saveRDS(Model_Correlation, MCMCglmm_filename)
i = i + 1
####################


#Univariate regression is covariance between response and predictor/variance in predictor
Model_Correlation_VCV<- mcmc.list(lapply(Model_Correlation, function(m) m$VCV))
plot(Model_Correlation_VCV) 
gelman.diag(Model_Correlation_VCV,multivariate = FALSE)

Model_Correlation_VCV_gg<- ggmcmc::ggs(Model_Correlation_VCV)

# for between species use the parameters ending in species 
Covariance<- Model_Correlation_VCV_gg %>%
  filter(Parameter == 'traitFissionBinary.2:traitEarlyGermlineBinary.2.species_rotl') %>% 
  rename(Covariance = value) %>% select(-Parameter)

#in previous model was: 'traitFissionBinary.1:traitEarlyGermlineBinary.1.species'

VariancePredictor<- Model_Correlation_VCV_gg %>%
  filter(Parameter == 'traitFissionBinary.2:traitFissionBinary.2.species_rotl') %>%
  rename(PredictorVariance = value) %>% select(-Parameter)

# in previous model was: traitFissionBinary.1:traitFissionBinary.1.species

Regression<- merge(Covariance, VariancePredictor) %>%
  mutate(value = Covariance / PredictorVariance)
glimpse(Regression)
ggplot(Regression, aes(x = Iteration, y = value, colour = as.factor(Chain))) + geom_line(size = 1, alpha = 0.5) + geom_smooth() + theme_minimal()
dev.off()
}


