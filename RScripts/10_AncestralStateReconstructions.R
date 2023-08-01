# Ancestral State Reconstructions using Charlie's Bayes functions
# load the files containing the functions
devtools::source_url("https://raw.githubusercontent.com/charliecornwallis/Rfunctions/master/ASR_functions.R")


source('RScripts/pMCMC_Script.R')

#tree<- read.tree("Data/Phylogeny_ROTL.txt") ## BEWARE HERE ABOUT THAT ANNOYING YAMADAELLA NAMe
#tree<- tr

#inv_tree<-inverseA(tree)$Ainv
#plot(tree)


# while testing:
# iterations<- 8000
# burnin<- 1000
# n_chains<- 6
# thinning <- 1

# Running models to reconstruct the traits

# cell types
#Read in the priors
dfACR<- df
#dfACR<- subset(dfACR, GermNumeric != 0)
#h dfACR<- subset(dfACR, !is.na(GermNumeric))



priors= list(p3,p2,p1)

i = 1
for(prior in priors){
  # run MCMCglmm
  print(prior)
  MCMCglmm_Types <- mclapply(1:n_chains, function(i){
    MCMCglmm(Types ~ 1, #-1 here removes the intercept equivalent to 0 in brms
             random = ~species_rotl, ginverse=list(species_rotl=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
             family ="poisson",data = dfACR,prior=prior, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=T, pl = T)
  }, mc.cores = n_cores)
  
  names(MCMCglmm_Types)<- c('chain1','chain2','chain3','chain4','chain5','chain6')

  MCMCglmm_filename<- paste(PathForAnalyses, 'R_Objects/ModelOutputs/p', i,  '/AncestralStates_CellTypes.RDS', sep = '')
  saveRDS(MCMCglmm_Types, MCMCglmm_filename)
  i = i + 1
  


#MCMCglmm_Types<- readRDS(MCMCglmm_filename)


predTypes<-pred_states_mcmcglmm(mr="No", # is the model multi-response?
                           trees=tr, # what is the phylogeny used?
                           phy_name="species_rotl", # the in-model name for the phylogeny?
                           model=MCMCglmm_Types$chain2, # the model name
                           dat=as.data.frame(dfACR), #the dataframe
                           trait1="Type", #the response variable
                           binomial="No") # whether or not this trait is binary -> since I’m reconstructing our continuous response variable, set to ‘no’
#state1="early", # name of first state of binary trait
#state2="late", # name of second state of trait
#cutoff=0.5)

print(predTypes)

}

# cell numbers
i = 1
for(prior_set in priors){
  # run MCMCglmm
  MCMCglmm_Number <- mclapply(1:n_chains, function(i){
    MCMCglmm(Number ~ 1, #-1 here removes the intercept equivalent to 0 in brms
             random = ~species_rotl, ginverse=list(species_rotl=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
             family ="poisson",data = dfACR,prior=prior_set, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=T, pl = T)
  }, mc.cores = n_cores)
  
  names(MCMCglmm_Number)<- c('chain1','chain2','chain3','chain4','chain5','chain6')
  
  MCMCglmm_filename<- paste(PathForAnalyses, 'R_Objects/ModelOutputs/p', i, '/AncestralStates_CellNumbers.RDS', sep = '')
  saveRDS(MCMCglmm_Number, MCMCglmm_filename)
  i = i + 1
  


#MCMCglmm_Number<- readRDS(MCMCglmm_filename)


predNumber<-pred_states_mcmcglmm(mr="No", # is the model multi-response?
                                trees=tr, # what is the phylogeny used?
                                phy_name="species_rotl", # the in-model name for the phylogeny?
                                model=MCMCglmm_Number$chain1, # the model name
                                dat=as.data.frame(dfACR), #the dataframe
                                trait1="Number", #the response variable
                                binomial="No") # whether or not this trait is binary -> since I’m reconstructing our continuous response variable, set to ‘no’
#state1="early", # name of first state of binary trait
#state2="late", # name of second state of trait
#cutoff=0.5)
predNumber$`log(estimate)`<- log(predNumber$estimate)
print(predNumber)
}


# fission
i = 1
  MCMCglmm_Fission <- mclapply(1:n_chains, function(i){
    MCMCglmm(as.factor(Fission) ~ 1, #-1 here removes the intercept equivalent to 0 in brms
             random = ~species_rotl, ginverse=list(species_rotl=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
             family ="categorical",data = dfACR,prior=p_binary, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=T, pl = T)
  }, mc.cores = n_cores)
  
  names(MCMCglmm_Fission)<- c('chain1','chain2','chain3','chain4','chain5','chain6')
  
  MCMCglmm_filename<- paste(PathForAnalyses, 'R_Objects/ModelOutputs/p5/AncestralStates_Fission.RDS', sep = '')
  saveRDS(MCMCglmm_Fission, MCMCglmm_filename)

  


#MCMCglmm_Fission<- readRDS(MCMCglmm_filename)


predFission<-pred_states_mcmcglmm(mr="No", # is the model multi-response?
                                 trees=tr, # what is the phylogeny used?
                                 phy_name="species_rotl", # the in-model name for the phylogeny?
                                 model=MCMCglmm_Fission$chain1, # the model name
                                 dat=as.data.frame(dfACR), #the dataframe
                                 trait1="Fission", #the response variable
                                 binomial="Yes", # whether or not this trait is binary -> since I’m reconstructing our continuous response variable, set to ‘no’
                                 state1="0", # name of first state of binary trait
                                 state2="1", # name of second state of trait
                                cutoff=0.5)

print(predFission)





dfACR$GermNumeric[dfACR$GermNumeric == 1] <- 'early'
dfACR$GermNumeric[dfACR$GermNumeric == 2] <- 'late'

# germlines
i = 1

  # run MCMCglmm
  MCMCglmm_Germline <- mclapply(1:n_chains, function(i){
    MCMCglmm(as.factor(GermNumeric) ~ 1, #-1 here removes the intercept equivalent to 0 in brms
             random = ~species_rotl, ginverse=list(species_rotl=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
             family ="categorical",data = dfACR,prior=p_binary, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=T, pl = T)
  }, mc.cores = n_cores)
  
  names(MCMCglmm_Germline)<- c('chain1','chain2','chain3','chain4','chain5','chain6')
  
  MCMCglmm_filename<- paste(PathForAnalyses, 'R_Objects/ModelOutputs/p5/AncestralStates_Germlines.RDS', sep = '')
  saveRDS(MCMCglmm_Germline, MCMCglmm_filename)
  i = i + 1
  
  
  
 # MCMCglmm_Germline<- readRDS(MCMCglmm_filename)
  
  predGermline<-pred_states_mcmcglmm(mr="No", # is the model multi-response?
                                    trees=tr, # what is the phylogeny used?
                                    phy_name="species_rotl", # the in-model name for the phylogeny?
                                    model=MCMCglmm_Germline$chain1, # the model name
                                    dat=as.data.frame(dfACR), #the dataframe
                                    trait1="GermNumeric", #the response variable
                                    binomial="Yes", # whether or not this trait is binary -> since I’m reconstructing our continuous response variable, set to ‘no’
                                    state1="early", # name of first state of binary trait
                                    state2="late", # name of second state of trait
                                    cutoff=0.5)
  
  print(predGermline)



colnames(predTypes)<- c('species','CellTypes')
colnames(predNumber)<- c('species','CellNumber','logNumber')
colnames(predGermline)<- c('species','GermlineEstimate','GermlineState')
colnames(predFission)<- c('species','FissionEstimate','FissionState')


preds<- cbind(predTypes, predNumber[,-1], predGermline[,-1], predFission[,-1])
preds$CellNumberBackTrans<- exp(preds$CellNumber)
preds$CellTypesBackTrans<- exp(preds$CellTypes)

preds$TypesLogNumber <- preds$CellTypesBackTrans / log(preds$CellNumberBackTrans)

predsNodes<- subset(preds, grepl('Node',preds$species))
predsNodes$node<- 1:tr$Nnode+Ntip(tr)

predsTips<- subset(preds, !grepl('Node',preds$species))






ph <- ggtree(tr) %<+% predsNodes
ph <- ph + geom_point(aes(size=TypesLogNumber, colour = GermlineState, shape = FissionState))
ph

gg_df<- dfACR %>%
  mutate(id = species_rotl) %>%
  select(c('id','GermNumeric','Types','Number','Fission'))
#colnames(gg_df)<- c('id','GermlineState','Types','Number','Fission')

ph2<- ph %<+% gg_df
ph2<- ph2 + geom_tippoint(aes(colour = GermNumeric, shape = Fission, size = Types/log(Number)))
ph2 <- ph2 + geom_tiplab()
ph2

ph2<- ph + ggtree::geom_facet(panel = 'Cell Types', data = gg_df, geom = geom_col, mapping = aes(x= Types, fill = (as.factor(GermNumeric))), orientation = 'y')
ph3<- ph2 + ggtree::geom_facet(panel = 'Log(Cell Number)', data = gg_df, geom = geom_col, mapping = aes(x= log(Number), fill = (as.factor(GermNumeric))), orientation = 'y')
ph4<- ph3 + ggtree::geom_facet(panel = 'Cell Types / log(Cell Number)', data = gg_df, geom = geom_col, mapping = aes(x= Types/log(Number), fill = (as.factor(GermNumeric))), orientation = 'y')

ph5<- ph4 %<+% gg_df + geom_tippoint(aes(colour = GermNumeric, shape = Fission, size = Types/log(Number)))
ph5

ph_simple<- ph + ggtree::geom_facet(panel = 'Cell Types / log(Cell Number)', data = gg_df, geom = geom_col, mapping = aes(x= (Types/log(Number)), fill = (as.factor(GermNumeric))), orientation = 'y')
ph_simple<- ph_simple %<+% gg_df + geom_tippoint(aes(colour = GermNumeric, shape = Fission, size = Types/log(Number)))

pdf(height=20, width=15, paste(PathForAnalyses, "Figures/SimplePhylogeny.pdf", sep = ''), useDingbats = F)
print(ph_simple)
dev.off()



ph <- ggtree(tr) %<+% predsNodes
ph <- ph + geom_point(aes(size=TypesLogNumber, colour = GermlineState, shape = FissionState))
ph


### drawing figure with cell numbers and types

ph2<- ph %<+% gg_df
ph2<- ph2 + geom_tippoint(aes(colour = GermNumeric, shape = Fission, size = Types/log(Number)))
#p2 <- p2 + geom_tiplab()

ph2<- ph + ggtree::geom_facet(panel = 'Cell Types', data = gg_df, geom = geom_col, mapping = aes(x= Types, fill = (as.factor(GermNumeric))), orientation = 'y')
ph3<- ph2 + ggtree::geom_facet(panel = 'Log(Cell Number)', data = gg_df, geom = geom_col, mapping = aes(x= log(Number), fill = (as.factor(GermNumeric))), orientation = 'y')
ph4<- ph3 + ggtree::geom_facet(panel = 'Cell Types / log(Cell Number)', data = gg_df, geom = geom_col, mapping = aes(x= Types/log(Number), fill = (as.factor(GermNumeric))), orientation = 'y')

ph5<- ph4 %<+% gg_df + geom_tippoint(aes(colour = GermNumeric, shape = Fission, size = Types/log(Number)))
ph5

ph6<- facet_widths(ph5, widths = c(1.5,1,1,1))

pdf(height=25, width=20, paste(PathForAnalyses, "Figures/AllPlotsPhylogeny.pdf", sep = ''), useDingbats = F)
print(ph6)
dev.off()
