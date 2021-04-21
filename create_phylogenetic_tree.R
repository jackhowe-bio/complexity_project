# script that creates phylogenetic tree file

#load packages
library(ape)
library(ggplot2)
library(brms)

# create tree


myTreeAnimals <- read.tree(text='
( 

  (
    (Hirudo medicinalis,
    Lumbricus terrestris,
    Nais variabilis,
    Dinophilus gyrociliatus,
    Diurodrilus westheidei,
    apodotrocha progenerans,
    Aeolosoma tenebrarum,
    Dasybranchus caducus,
    Pisione remota,
    Pomatoceros triqueter),
  
    (Periplaneta americana,
    Callinectes sapidus),
  )
  
  pycnophyes frequens,

  (Morone saxatilis,
  Salmo gairdneri,
  Canis familiaris,
  Mus musculus),

  (Hydra vulgaris,
  microhydra rideri,
  cyanea cyanea,
  haliclystus haliclystus),

  (Pleurobrachia),

  (conocyema polymorpha,
  dicyema typhus,
  dicyemmenea abelis,
  dicyemmenea lameerei),

  (turbanella cornuta,
  chordodasys antennatus),

  (valvognathia pogonostoa,
  rastrognathhia macrostoma),

  (loxosoma sultana,
  pedicillina echinata),

  (amphibola crenata,
  neomenia carinata),

  (caenorhabditis elegans,
  rhabditis monhystera,)

  (rhopalura granosa),
  
  (Trichoplax adhaerens),
  
  macrostomum gigas,
  entersostomula graffi,
  dugesia mediterranea,
  spongilla lacustris,
  notholca acuminata,
  apsilus vorax,
  anaperus sulcatus
  

  
);')
plot(myTreeAnimals)


library(rotl)
taxa <- c("Hyla", "Salmo", "Diadema", "Nautilus")
resolved_names <- tnrs_match_names(taxa)

df<- read.csv('germline_data_1.0.csv')
animals <- subset(df, kingdom == 'Animalia')
resolved_names_animals <- tnrs_match_names(animals$Species, context_name = 'Animals')


