# script that creates phylogenetic tree file

#load packages
library(ape)
library(ggplot2)
library(brms)
library(rotl) #see https://cran.r-project.org/web/packages/rotl/vignettes/rotl.html


df<- read.csv('germline_data_1.0.csv')
df[nrow(df),]
df<- df[-197,]
df[nrow(df),]

# create tree
#ALL 
ResolvedNames <- tnrs_match_names(df$species.updated.rotl, context_name = 'All life')
ResolvedNames$IsInTree <- is_in_tree(ResolvedNames$ott_id)
ResolvedNamesInTree<- subset(ResolvedNames, IsInTree==T)

AllTree<- tol_induced_subtree(ResolvedNamesInTree$ott_id)

# draw tree
plot(AllTree, no.margin = TRUE, cex = 0.5, label.offset = 0.5)

#save tree
write.tree(AllTree, file='phylogeny_all.txt')

#resolve polytomies
ResolvedPolytomiesTree<- multi2di(AllTree)
plot(ResolvedPolytomiesTree, no.margin = TRUE, cex = 0.5, label.offset = 0.5)
write.tree(AllTree, file='phylogeny_all_res_polytomy.txt')


write.csv(ResolvedNames, 'phylogeny_species_names.csv')
write.csv(ResolvedNamesInTree, 'phylogeny_species_names_in_tree.csv')
