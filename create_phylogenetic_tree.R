# script that creates phylogenetic tree file

#load packages
library(ape)
library(ggplot2)
library(brms)
library(rotl) #see https://cran.r-project.org/web/packages/rotl/vignettes/rotl.html


df<- read.csv('germline_data_1.0.csv')
df<- unique(df)
# create tree

#ALL 
ResolvedNames <- tnrs_match_names(df$species.updated.rotl, context_name = 'All life')
ResolvedNames <- subset(ResolvedNames, !is.na(ott_id))

ResolvedNames$IsInTree <- is_in_tree(ResolvedNames$ott_id)
ResolvedNamesInTree<- subset(ResolvedNames, IsInTree==T)
AllTree<- tol_induced_subtree(ResolvedNamesInTree$ott_id)
plot(AllTree, no.margin = TRUE)


