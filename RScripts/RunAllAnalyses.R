# data import using rotl project
setwd('/Users/pcx971/Documents/oxford/complexity/complexity_project/')# See https://cran.r-project.org/web/packages/rotl/vignettes/meta-analysis.html 
# Link provides a vignette for running comparative analysis using rtol project

# load the library
library(rotl)
library(ape)
library(tidyverse)
library(MCMCglmm)
library(parallel)
library(coda)
library(ggthemes)
library(treeio)
library(ggpubr)
library(bayesplot)
library(ggtree)
library(ggplot2)
library(binom)
library(ggbeeswarm)

# make the structure for the directory
system("mkdir -p Results/{SpeciesLevel,GenusLevel}/{All,OnlyAnimals,NoAnimals}/R_Objects/Model{Outputs,Diagnostics}/{p1,p2,p3,p4,p5}")
system("mkdir -p Results/{SpeciesLevel,GenusLevel}/{All,OnlyAnimals,NoAnimals}/Figures")
system("mkdir -p Results/{SpeciesLevel,GenusLevel}/{All,OnlyAnimals,NoAnimals}/Data/Anvio")

#import the data and get it into the right format
dataframe_input<- 
  read_csv('Data/germline_data_1.5.csv', col_types = list( #import data with defining column types
    cell_types = col_number(),
    cell_number = col_number(),
    clonal = col_factor(),
    SexualReproductionObserved_Species = col_factor(),
    AsexualGameticObserved_Species = col_factor(),
    FissionOrBuddingObserved_Species = col_factor(),
    SexualReproductionObserved_Genus = col_factor(),
    AsexualGameticObserved_Genus = col_factor(),
    FissionOrBuddingObserved_Genus = col_factor(),
    #fragmentation_or_budding = col_factor(),
    germline_timing = col_factor()
  ))


# set the priors to be used
p1=list(R = list(V = 1, nu=0.002), 
        G = list(G1=list(V = 1, nu = 0.002)))

p2=list(R=list(V=1, nu=1), 
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.var=1000)))

p3=list(R=list(V=1, nu=2),
        G=list(G1=list(V=1, nu=2, alpha.mu=0, alpha.var=1000)))

p4=list(B=list(mu=c(0,0), V=diag(c(1+pi^2/3,1+pi^2/3))),
        R = list(V = diag(2),nu=1, fix=1), #Residual variance not identifiable for binary variables
        G = list(G1=list(V = diag(2), nu = 1, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000)))))#parameter expanded priors, usually good for binary data

# used for ancestral state reconstructions
p_binary =
  list(B=list(mu=c(0), V=1+pi^2/3),R = list(V = 1, nu=1,fix=1), G = list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))

# set other MCMCglmm parameters
iterations<- 8000000
burnin<- 1000000
thinning <- 1000
n_chains<- 6

n_cores<- 4

# # while testing:
# iterations<- 8000
# burnin<- 1000
# n_chains<- 6
# thinning <- 1

# objects to keep = 
objects_to_keep = ls()
objects_to_keep[length(objects_to_keep) + 1]<- "objects_to_keep"

############################################################
#Species Level
############################################################
  
      #########################
      # All
      #########################
 PathForAnalyses<- "Results/SpeciesLevel/All/"
#
# # this is where it needs subset
 multicell_data<- dataframe_input %>%
   subset(germline_timing != '0' & germline_timing != '?')
#
 # then prepare the dataset and the tree: this saves to the file that is loaded next
 source("RScripts/1_DataPrep.R")
#
 #Read in the data and the tree
 File=paste(PathForAnalyses, 'R_Objects/metadata.RDS', sep = '')
 df = readRDS(File)
 File=paste(PathForAnalyses, 'R_Objects/inv_tree.RDS', sep = '')
 inv_tree = readRDS(File)
#
 #the optimisation needs run once
 source("RScripts/3_Optimisation.R")

# # run the models
#
 source("RScripts/4_Model1_ROTL.R")
 source("RScripts/5_Model2_ROTL.R")
 source("RScripts/6_Model3_ROTL.R")
 source('RScripts/7_Model4_ROTL.R')
#
 # then the correlation models
 source('RScripts/8_Multiresponse_ROTL.R')
 source('RScripts/9_GermlineWithoutCellNumber.R')
 source('RScripts/10_AncestralStateReconstructions.R')
#
# # draw the figures
 source('RScripts/11_Figures.R')
#
 objects_to_remove = ls()[(ls() %in% objects_to_keep) == F]
 rm(list = objects_to_remove)
#
#
#       ################################
#       # Without Animals
#       ################################
#
 PathForAnalyses<- "Results/SpeciesLevel/NoAnimals/"
#
 # this is where it needs subset
 multicell_data<- dataframe_input %>%
   subset(germline_timing != '0' & germline_timing != '?') %>%
   subset(kingdom != 'Animalia')

 # then prepare the dataset and the tree: this saves to the file that is loaded next
 source("RScripts/1_DataPrep.R")

 #Read in the data and the tree
 File=paste(PathForAnalyses, 'R_Objects/metadata.RDS', sep = '')
 df = readRDS(File)
 File=paste(PathForAnalyses, 'R_Objects/inv_tree.RDS', sep = '')
 inv_tree = readRDS(File)

 # run the models
 source("RScripts/3_Optimisation.R")
 source("RScripts/4_Model1_ROTL.R")
 source("RScripts/5_Model2_ROTL.R")
 source("RScripts/6_Model3_ROTL.R")
 source('RScripts/7_Model4_ROTL.R')
#
# # then the correlation models
 source('RScripts/8_Multiresponse_ROTL.R')
 source('RScripts/9_GermlineWithoutCellNumber.R')
#
 # and ancestral state reconstructions
 source('RScripts/10_AncestralStateReconstructions.R')
#
# # draw the figures
 source('RScripts/11_Figures.R')
#
# # then tidy up the objects in memory
 objects_to_remove = ls()[(ls() %in% objects_to_keep) == F]
 rm(list = objects_to_remove)
#
#
#     ##################################
#     # Only animals
#     ##################################
#
 PathForAnalyses<- "Results/SpeciesLevel/OnlyAnimals/"
#
 # this is where it needs subset
 multicell_data<- dataframe_input %>%
   subset(germline_timing != '0' & germline_timing != '?') %>%
   subset(kingdom == 'Animalia')

# # then prepare the dataset and the tree: this saves to the file that is loaded next
 source("RScripts/1_DataPrep.R")

# #Read in the data and the tree
 File=paste(PathForAnalyses, 'R_Objects/metadata.RDS', sep = '')
 df = readRDS(File)
 File=paste(PathForAnalyses, 'R_Objects/inv_tree.RDS', sep = '')
 inv_tree = readRDS(File)

 # run the models
 source("RScripts/3_Optimisation.R")
 source("RScripts/4_Model1_ROTL.R")
 source("RScripts/5_Model2_ROTL.R")
 source("RScripts/6_Model3_ROTL.R")
 source('RScripts/7_Model4_ROTL.R')

 # then the correlation models
 source('RScripts/8_Multiresponse_ROTL.R')
 source('RScripts/9_GermlineWithoutCellNumber.R')

 # and ancestral state reconstructions
 source('RScripts/10_AncestralStateReconstructions.R')

 # draw the figures
 source('RScripts/11_Figures.R')

 # then tidy up the objects in memory
 objects_to_remove = ls()[(ls() %in% objects_to_keep) == F]
 rm(list = objects_to_remove)



##########################################
# Genus Level
#########################################
    ##########################
    ##### All
    #########################
PathForAnalyses<- "Results/GenusLevel/All/"

# this is where it needs subset
multicell_data<- dataframe_input %>%
  subset(germline_timing != '0' & germline_timing != '?')

# then prepare the dataset and the tree: this saves to the file that is loaded next
source("RScripts/1_DataPrep.R")

#Read in the data and the tree
File=paste(PathForAnalyses, 'R_Objects/metadata.RDS', sep = '')
df = readRDS(File)
File=paste(PathForAnalyses, 'R_Objects/inv_tree.RDS', sep = '')
inv_tree = readRDS(File)

# # run the models
# #source("RScripts/3_Optimisation.R")
# source("RScripts/4_Model1_ROTL.R")
# source("RScripts/5_Model2_ROTL.R")
# source("RScripts/6_Model3_ROTL.R")
# source('RScripts/7_Model4_ROTL.R')
# 
# # then the correlation models
# source('RScripts/8_Multiresponse_ROTL.R')
# source('RScripts/9_GermlineWithoutCellNumber.R')
# 
# # and ancestral state reconstructions
# source('RScripts/10_AncestralStateReconstructions.R')

# draw the figures
source('RScripts/11_Figures.R')

# then tidy up the objects in memory
objects_to_remove = ls()[(ls() %in% objects_to_keep) == F]
rm(list = objects_to_remove)




################################
# Without Animals
################################

PathForAnalyses<- "Results/GenusLevel/NoAnimals/"

# this is where it needs subset
multicell_data<- dataframe_input %>%
  subset(germline_timing != '0' & germline_timing != '?') %>%
  subset(kingdom != 'Animalia')

# then prepare the dataset and the tree: this saves to the file that is loaded next
source("RScripts/1_DataPrep.R")

#Read in the data and the tree
File=paste(PathForAnalyses, 'R_Objects/metadata.RDS', sep = '')
df = readRDS(File)
File=paste(PathForAnalyses, 'R_Objects/inv_tree.RDS', sep = '')
inv_tree = readRDS(File)

# run the models
#source("RScripts/3_Optimisation.R")
source("RScripts/4_Model1_ROTL.R")
source("RScripts/5_Model2_ROTL.R")
source("RScripts/6_Model3_ROTL.R")
source('RScripts/7_Model4_ROTL.R')

# then the correlation models
source('RScripts/8_Multiresponse_ROTL.R')
source('RScripts/9_GermlineWithoutCellNumber.R')

# and ancestral state reconstructions
source('RScripts/10_AncestralStateReconstructions.R')

# draw the figures
source('RScripts/11_Figures.R')
# then tidy up the objects in memory
objects_to_remove = ls()[(ls() %in% objects_to_keep) == F]
rm(list = objects_to_remove)



##################################
# Only animals
##################################

PathForAnalyses<- "Results/GenusLevel/OnlyAnimals/"

# this is where it needs subset
multicell_data<- dataframe_input %>%
  subset(germline_timing != '0' & germline_timing != '?') %>%
  subset(kingdom == 'Animalia')

# then prepare the dataset and the tree: this saves to the file that is loaded next
source("RScripts/1_DataPrep.R")

#Read in the data and the tree
File=paste(PathForAnalyses, 'R_Objects/metadata.RDS', sep = '')
df = readRDS(File)
File=paste(PathForAnalyses, 'R_Objects/inv_tree.RDS', sep = '')
inv_tree = readRDS(File)

# run the models
#source("RScripts/3_Optimisation.R")
source("RScripts/4_Model1_ROTL.R")
source("RScripts/5_Model2_ROTL.R")
source("RScripts/6_Model3_ROTL.R")
source('RScripts/7_Model4_ROTL.R')

# then the correlation models
source('RScripts/8_Multiresponse_ROTL.R')
source('RScripts/9_GermlineWithoutCellNumber.R')

# and ancestral state reconstructions
source('RScripts/10_AncestralStateReconstructions.R')

# draw the figures
source('RScripts/11_Figures.R')

# then tidy up the objects in memory
objects_to_remove = ls()[(ls() %in% objects_to_keep) == F]
rm(list = objects_to_remove)



