# script that creates phylogenetic tree file

#load packages
library(ape)
library(ggplot2)
library(brms)
library(rotl) #see https://cran.r-project.org/web/packages/rotl/vignettes/rotl.html


df<- read.csv('data/germline_data_1.1.csv')
df[nrow(df),]
df<- df[-197,]
df[nrow(df),]

# create tree
#ALL 
ResolvedNames <- tnrs_match_names(df$species.updated.rotl, context_name = 'All life')
ResolvedNames$IsInTree <- is_in_tree(ResolvedNames$ott_id)
ResolvedNamesInTree<- subset(ResolvedNames, IsInTree==T)

AllTree<- tol_induced_subtree(ResolvedNamesInTree$ott_id, label_format = 'id')


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

#CKC edits
#Can also produce a tree using taxonomic classifications - here's an example

species<- df %>% select(order,family,genus,species) %>% arrange(order,family,genus,species)

#Add outgroup
out<-data.frame(order="Ulvales",family="Kornmanniaceae",genus="Pseudendoclonium",species="Pseudendoclonium_basiliense")
species<-rbind(species,out)
species<- species %>% mutate(across(,~as.factor(.x)))

tree<-as.phylo(~order/family/genus/species, data=species)
tree<-multi2di(tree, random=T)
tree<-makeNodeLabel(tree, method="number", prefix="Node")
tree<-root(tree, outgroup="Pseudendoclonium_basiliense", resolve.root=T)
tree<-compute.brlen(tree, method = "Grafen", power = 1)
tree<-di2multi(tree, tol=1e-25)
tree<-chronoMPL(tree)

plot(tree)
is.binary.tree(tree)
is.ultrametric(tree)


