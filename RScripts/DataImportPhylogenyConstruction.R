# script that imports the data files, constructs the phylogenies,
#and the inverse matrix for MCMCglmm functions

library(ape)
library(tidyverse)
library(MCMCglmm)
library(parallel)
library(coda)
library(ggthemes)

df<- read_csv('data/germline_data_1.2.csv', col_types = list( #import data with defining column types
  cell_types = col_integer(),
  cell_number = col_integer(),
  clonal = col_factor(),
  SexualReproductionObserved_Species = col_factor(),
  AsexualGameticObserved_Species = col_factor(),
  FissionOrBuddingObserved_Species = col_factor(),
  SexualReproductionObserved_Genus = col_factor(),
  AsexualGameticObserved_Genus = col_factor(),
  FissionOrBuddingObserved_Genus = col_factor(),
  fragmentation_or_budding = col_factor(),
  germline_timing = col_factor()
))

#create a genus column
df<- df %>% 
  separate(species, into = c('genus', 'species_name'), remove = F) 

df$cell_number<- as.integer(as.numeric(df$cell_number))

early = c('1','1,2','2')

df$germline_timing_simple<- ifelse(df$germline_timing %in% early, 'early', df$germline_timing)  # early germline segregation
df$germline_timing_simple<- ifelse(df$germline_timing == '0', 'no_germline', df$germline_timing_simple) # no germline/soma divide
df$germline_timing_simple<- ifelse(df$germline_timing == '3', 'adult', df$germline_timing_simple) # germline segregated in adult tissues

df<- subset(df, !is.na(cell_number))

#Get species in dataset and order by classification
species<-df %>% select(kingdom, phylum, class, order,taxon,genus,species) %>% arrange(kingdom, phylum, class, order,taxon,genus,species)

species<- species %>% mutate(across(.cols = everything(),~as.factor(.x))) #make factors
tree<-as.phylo(~kingdom/phylum/class/order/taxon/genus/species, data=species) #draw tree
tree<-multi2di(tree, random=T) #collapse and resolve multichotomies (randomly)
tree<-makeNodeLabel(tree, method="number", prefix="Node") #make node labels
tree<-compute.brlen(tree, method = "Grafen", power = 1) #compute branch lengths
tree<-di2multi(tree, tol=1e-25) # collapse and resolve multichotomies (again?)
tree<-chronoMPL(tree) #estimates age of nodes from distance to all tips
treeTrimmed<- drop.tip(tree, which(is.na(tree$tip.label))) #drop any tips where the label is NA (not sure this is an issue any longer)

#should we include an outgroup? 

is.binary(treeTrimmed)
is.ultrametric(treeTrimmed)

pdf('figures/phylogenies.pdf')
plot(treeTrimmed, cex = 0.4)
dev.off()

#Create inverse for running MCMCglmm models 
inv_tree<-inverseA(treeTrimmed)$Ainv

df <- df %>%
  filter(df$species %in% dimnames(inv_tree)[[1]])
df<- as.data.frame(df)