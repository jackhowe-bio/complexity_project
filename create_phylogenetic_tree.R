# script that creates phylogenetic tree file

#load packages
library(ape)
library(ggplot2)
library(brms)
library(rotl) #see https://cran.r-project.org/web/packages/rotl/vignettes/rotl.html


df<- read.csv('germline_data_1.0.csv')

# create tree

#animals
animals <- subset(df, kingdom == 'Animalia')
ResolvedNamesAnimals <- tnrs_match_names(animals$species.updated, context_name = 'Animals')

AnimalTree<- tol_induced_subtree(ResolvedNamesAnimals$ott_id[-c(33,41)])
plot(AnimalTree, no.margin = TRUE)
