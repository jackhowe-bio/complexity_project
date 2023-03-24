# interpreting model 1 and 2to get statistics

setwd('/Users/pcx971/Documents/oxford/complexity/complexity_project/HPC_Analyses/')

## load packages
library(ape)
library(tidyverse)
library(MCMCglmm)
library(parallel)
library(coda)
library(ggthemes)
library(ggpubr)


for(Model in c(1,2)){
  for(prior_set in c('p1','p2','p3')){
    
    #read in model
    input_filename<- paste('RScripts/ModelOutputs/',prior_set,'/Model',Model,'_ROTL_',prior_set,'.RDS', sep = '')
    ModelOutput <- readRDS(input_filename)
  
    Sol_mcmc<- mcmc.list(lapply(ModelOutput, function(m) m$Sol))
    Sol_gg<- ggmcmc::ggs(Sol_mcmc)  
    Sol_ggSubset<- subset(Sol_gg,  grepl('Fission', Sol_gg$Parameter)) #need to subset so we don't get a plot for each species
    
    Sol_traceplot <- ggmcmc::ggs_traceplot(Sol_ggSubset) + theme_minimal()
    Sol_caterpillar<- ggmcmc::ggs_caterpillar(Sol_ggSubset) + theme_minimal()
    
    SolGelman<- gelman.diag(Sol_mcmc,multivariate = FALSE)
    
    NoFission<- Sol_gg %>%
      filter(Parameter == 'Fission0') %>% 
      rename(Fission0 = value) %>% select(-Parameter)
    
    Fission<- Sol_gg %>%
      filter(Parameter == 'Fission1') %>%
      rename(Fission1 = value) %>% select(-Parameter)
    
    Difference<- merge(Fission, NoFission) %>%
      mutate(Difference = Fission1 - Fission0)
  
    HDI_90<- bayestestR::ci(Difference$Difference, method = 'HDI', ci = 0.90)
    HDI_95<- bayestestR::ci(Difference$Difference, method = 'HDI', ci = 0.95)
    
    DifPlotM<- ggplot(Difference, aes(colour = as.factor(Chain), x = Difference)) + geom_density() + theme_minimal() + geom_segment(x = HDI_90$CI_low, y = 0, xend = HDI_90$CI_high, yend = 0, colour = 'grey', size = 4) + geom_segment(x = HDI_95$CI_low, y = 0, xend = HDI_95$CI_high, yend = 0, colour = 'grey', size = 1)
    
    SummaryPlot<- ggarrange(Sol_traceplot, Sol_caterpillar, DifPlotM, labels = "AUTO")
    SummaryPlot
    
    
    OutputFigureName<- paste('RScripts/ModelOutputs/',prior_set,'/Model',Model,'SummaryFigure_',prior_set,'.pdf',sep = '')
    pdf(OutputFigureName)
    print(SummaryPlot)
    dev.off()
    }  
}  
