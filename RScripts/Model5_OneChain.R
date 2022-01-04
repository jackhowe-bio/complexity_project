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


p1=list(R = list(V = 1, nu=0.002)) #sets prior for residual variance, the defaults are used as priors for fixed effects (see MCMCglmm course notes)

M1 <-   MCMCglmm(cell_number ~ FissionOrBuddingObserved_Species-1, #-1 here removes the intercept equivalent to 0 in brms
           family ="poisson",data = df,prior=p1, nitt=300000, burnin=3000, thin=thinning,verbose = T,pr=T)



####SAVE OUTPUT #### 
#save image of environment so it doesn't need run each time
#save.image(file = 'RScripts/Model5.Rdata')

#CKC Processing MCMCglmm models ####

pacman::p_load(devtools)
#CKC - function for processing MCMCglmm models
source('RScripts/pMCMCglmmScript.R')

#Model formula being processed
#MCMCglmm(cell_number ~ FissionOrBuddingObserved_Species-1...

#Extract 1 model
M1<-M1

#Apply function (description of function terms can be found on github: https://github.com/charliecornwallis/Rfunctions/blob/master/MCMCglmmProc.R)
SItablesXL<-MCMCglmmProc(model=M1,pvalues="all",link="poisson",start_row=1,workbook=SItablesXL,create_sheet="yes",sheet="Table S1",title="Table S1: Jack project", fixed_names=c("FissionOrBuddingObserved_Species0","FissionOrBuddingObserved_Species1","FissionOrBuddingObserved_Species?","FissionOrBuddingObserved_Species"),
                         fixed_diffinc = c("all"),
                         Include_random = "yes",
                         randomvar_names=c("Residual"),
                         variances=c("units"),padding=3)

#The output is in excel format and can be written to excel using the openxlsx package

#Can also be converted to Rmd friendly format
#function for extracting df from xl workbook
xl_2_df = function(xltab,sheet=NULL){
  df<-readWorkbook(xltab,sheet=sheet)
  colnames(df)<-df[1,]
  df<-df %>% filter(pMCMC != "" & row_number() != 1)
  rownames(df)<-NULL
  return(df)
}

S1<-xl_2_df(xltab=SItablesXL,sheet="Table S1")
md_table(S1)
