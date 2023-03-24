# data import using rotl project
setwd('/Users/pcx971/Documents/oxford/complexity/complexity_project/HPC_Analyses/')# See https://cran.r-project.org/web/packages/rotl/vignettes/meta-analysis.html 
# Link provides a vignette for running comparative analysis using rtol project

# load the library
library(rotl)
library(ape)
library(tidyverse)
library(MCMCglmm)
library(parallel)
library(coda)
library(ggthemes)
#library(treeio)


#import the data and get it into the right format
multicell_data<- 
  read_csv('Data/germline_data_1.4.csv', col_types = list( #import data with defining column types
  cell_types = col_number(),
  cell_number = col_number(),
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

# take only the obligate
multicell_data<- subset(multicell_data, Obligate_or_facultative == 'obligate')

# generate genus and species names
multicell_data<- multicell_data %>% 
  separate(species, into = c('genus', 'species_name'), remove = F) 

#round cell number to whole number
multicell_data$cell_number<- ceiling(multicell_data$cell_number)

# change germline timing to 'early','adult','no_germline'
early = c('1','1,2','2')
multicell_data$germline_timing_simple<- ifelse(multicell_data$germline_timing %in% early, 'early', NA)  # early germline segregation
multicell_data$germline_timing_simple<- ifelse(multicell_data$germline_timing == '0', 'no_germline', multicell_data$germline_timing_simple) # no germline/soma divide
multicell_data$germline_timing_simple<- ifelse(multicell_data$germline_timing == '3', 'adult', multicell_data$germline_timing_simple) # germline segregated in adult tissues

#take out the NAs
multicell_data<- subset(multicell_data, !is.na(cell_number))
multicell_data<- subset(multicell_data, !is.na(cell_types))


# tidy up 
metadata<- multicell_data %>% select(species, kingdom, phylum, Obligate_or_facultative, Simple_or_complex, Sterile_soma,
                                           Proportion_of_sterile_cells, cell_types, cell_number, SexualReproductionObserved_Species,
                                           FissionOrBuddingObserved_Species, fragmentation_or_budding, germline_timing, germline_timing_simple) %>%
  mutate(cell_type_over_number = cell_types/ cell_number, somatic_cells = cell_types - 1, somatic_cells_over_number = somatic_cells/cell_number) %>%
  mutate(germline_timing_numeric = case_when(
    germline_timing_simple == 'early' ~ "1",
    germline_timing_simple == 'adult' ~ '2',
    germline_timing_simple == 'no_germline' ~ '0',
    germline_timing_simple == 3 ~ 'NA') ) 

colnames(metadata) <- c("species_name","Kingdom","Phylum","Oblig/Fac","Simple/Complex","SterileCells", "SterileProp","Types","Number","Sex",
                        "Fission","Fragmentation","GermTime","GermTimeSimp",
                        "Type/Number","SomaTypes","Soma/Number","GermNumeric")

metadata <- metadata %>%
  filter(!is.na(Fission))

write.csv(metadata, 'Data/AllMetadata.txt')

## GENERATE THE PHYLOGENY
#find the names from rotl
taxa <- tnrs_match_names(unique(metadata$species_name))
taxon_map <- structure(taxa$search_string, names = taxa$unique_name)
metadata$species_rotl<- taxa$search_string

# generate a tree
tr <- tol_induced_subtree(ott_id(taxa)[is_in_tree(ott_id(taxa))])

# make the tip labels readable for next steps
otl_tips <- strip_ott_ids(tr$tip.label, remove_underscores = TRUE)
tr$tip.label <- taxon_map[ otl_tips ]
names(tr$tip.label)<- NULL
#tr$node.label <- strip_ott_ids(tr$node.label, remove_underscores = TRUE) #not sure if necessary
original_tr<- tr
# remove the node labels (not necessary?)
#tr$node.label <- NULL


# collapse out the multichotomies, calculate branch lengths
tr<-multi2di(tr, random=T) #collapse and resolve multichotomies (randomly)
#tree<-makeNodeLabel(tree, method="number", prefix="Node") #make node labels
tr<-compute.brlen(tr, method = "Grafen", power = 1) #compute branch lengths
tr<-di2multi(tr, tol=1e-25) # collapse and resolve multichotomies (again?)
tr<-chronoMPL(tr) #estimates age of nodes from distance to all tips
tr<- drop.tip(tr, which(is.na(tr$tip.label))) #drop any tips where the label is NA (not sure this is an issue any longer)
AnvioTree<- tr
tr$node.label<- NULL 

#Create inverse for running MCMCglmm models 
inv_tree<-inverseA(tr)$Ainv

pdf('figures/phylogeniesROTL.pdf')
plot(tr, cex = 0.2)
dev.off()


# only take the rows that are in the tree
metadata <- metadata[metadata$species_rotl %in% tr$tip.label, ]

## metadata files
metadata$species_simple_underscores <- gsub(' ','_', metadata$species_rotl)
metadata<- as.data.frame(metadata)

#multicell_data<- subset(multicell_data, !is.na())

# write table
write.table(metadata, 'data/ROTL_metadata.txt', sep = '\t', quote = F, row.names = F, fileEncoding="UTF-8")
write.tree(AnvioTree, file = 'anvio/Phylogeny_ROTL.txt')


# save the inverse tree object for loading in future
saveRDS(inv_tree, file = 'RScripts/R_Objects/inv_tree.RDS')
saveRDS(metadata, file = 'RScripts/R_Objects/metadata.RDS')



anvio_file<- metadata %>%
  select(species_simple_underscores, species_rotl, Fission, GermTimeSimp, Number, Types) %>%
  mutate(TypesControlledForNumber = Types/log(Number))

write_tsv(anvio_file, "anvio/metadata.txt")
