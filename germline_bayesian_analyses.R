## ----packages---------------------------------------------------------------------------------------------------------------------------------------------------
library(ape)
library(ggplot2)
library(tidyverse)
library(knitr)
library(brms) #https://rdrr.io/cran/brms/f/vignettes/brms_phylogenetics.Rmd
library(tidybayes)
library(rotl) #see https://cran.r-project.org/web/packages/rotl/vignettes/rotl.html
library(phytools)
library(coda)
library(MCMCglmm)

knitr::purl('germline_bayesian_analyses.Rmd','germline_bayesian_analyses.R') #read this rmd to an .R file


## ----data_input-------------------------------------------------------------------------------------------------------------------------------------------------
df<- read.csv('data/germline_data_1.2.csv')
#create a genus column
df<- df %>% 
  separate(species.updated.rotl, into = c('genus', 'species'), remove = F) 
glimpse(df)

df$cell_number<- as.integer(as.numeric(df$cell_number))


## ----phylogeny_creation, message=F, results='hide'--------------------------------------------------------------------------------------------------------------
ResolvedNames <- tnrs_match_names(df$species.updated.rotl, context_name = 'All life') #search for similar names in the 'open tree of life' project (ROTL)cc
ResolvedNames$IsInTree <- is_in_tree(ResolvedNames$ott_id) #T/F, did the above find a match that can be put in the phylogeny? 
ResolvedNamesInTree<- subset(ResolvedNames, IsInTree==T) #subset to only those where present in the phylogeny

AllTree<- tol_induced_subtree(ResolvedNamesInTree$ott_id, label_format = 'id') #draw the phylogeny with ids
AllTreeNames<- tol_induced_subtree(ResolvedNamesInTree$ott_id, label_format = 'name') #draw the phylogeny with names

# tree with resolved polytomies:
ResolvedPolytomiesTree<- multi2di(AllTree, random=T)
ResolvedPolytomiesTreeNames<- multi2di(AllTreeNames, random=T)

plot(ResolvedPolytomiesTreeNames)

# write to files
write.tree(ResolvedPolytomiesTree, file='data/phylogeny_all.txt') #phylogeny
write.tree(ResolvedPolytomiesTree, file='data/phylogeny_all_res_polytomy.txt') #phylogeny with resolved polytomies
write.csv(ResolvedNames, 'data/phylogeny_species_names.csv') #list of names
write.csv(ResolvedNamesInTree, 'data/phylogeny_species_names_in_tree.csv') #list of names that are present in the tree


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get species in dataset and order by calssifications
species<-df %>% select(kingdom, phylum, class, order,family,genus,species.updated.rotl) %>% arrange(kingdom, phylum, class, order,family,genus,species.updated.rotl)

#Add outgroup
out<-data.frame(kingdom = 'Plantae', phylum = 'Chlorophyta', class = 'Ulvophyceae', order="Ulvales",family="Kornmanniaceae",genus="Pseudendoclonium",species.updated.rotl="Pseudendoclonium basiliense")
species<-rbind(species,out)
species<- species %>% mutate(across(,~as.factor(.x)))

tree<-as.phylo(~kingdom/phylum/class/order/family/genus/species.updated.rotl, data=species)
tree<-multi2di(tree, random=T)
tree<-makeNodeLabel(tree, method="number", prefix="Node")
#tree<-root(tree, outgroup="Pseudendoclonium basiliense", resolve.root=T)
tree<-compute.brlen(tree, method = "Grafen", power = 1)
tree<-di2multi(tree, tol=1e-25)
tree<-chronoMPL(tree)
treeTrimmed<- drop.tip(tree, which(is.na(tree$tip.label)))

plot(treeTrimmed)
is.binary.tree(treeTrimmed)
is.ultrametric(treeTrimmed)

#need to make tip labels match to be able to compare...
treeTrimmed$tip.label
ResolvedPolytomiesTreeNames$tip.label <- gsub("_"," ",ResolvedPolytomiesTreeNames$tip.label)
ResolvedPolytomiesTreeNames$tip.label <- str_remove(ResolvedPolytomiesTreeNames$tip.label, " \\(.*\\)$")

#Can then compare the trees using a cophylogeny plot
Cophy<-cophylo(ResolvedPolytomiesTreeNames,treeTrimmed,rotate=T)
pdf('figures/phylogenies.pdf')
plot.phylo(ResolvedPolytomiesTree, cex = 0.5)
plot(treeTrimmed, cex = 0.4)
dev.off()

pdf('figures/cophylogeny.pdf')
plot(Cophy,link.type="curved",link.lwd=3,link.lty="solid",
     link.col=make.transparent("blue",0.25),fsize=0.2)
dev.off()


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------
#read in phylogeny
phylo<- ape::read.tree('data/phylogeny_all.txt')
#subset df above with just those that are in the tree (note some missing)
names_in_tree<- read.csv('data/phylogeny_species_names_in_tree.csv')
names_in_tree$ott_id<- paste('ott',names_in_tree$ott_id, sep = '')
#phylo$tip.label %in% names_in_tree$ott_id #to check whether they're all present
names_in_tree<- subset(names_in_tree, names_in_tree$ott_id %in% phylo$tip.label)
#which species in the table are in the tree
df$species.updated.rotl<- tolower(df$species.updated.rotl)
df<- subset(df, tolower(df$species.updated.rotl) %in% tolower(names_in_tree$search_string))
#give them the ids in a column sot hat brms can match it up
df$species_id<- names_in_tree$ott_id


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# set branch lengths to 1 for covariance matrix
phylo_1b <- compute.brlen(phylo, 1)
#create covariance matrix
CovarMatrix <- ape::vcv.phylo(phylo_1b)


## ----plot_phylogeny---------------------------------------------------------------------------------------------------------------------------------------------
#plot phylogeny to check
plot(AllTree, no.margin = TRUE, cex = 0.5, label.offset = 0.5)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------
df$FissionOrBuddingObserved_Genus_nominal <- ifelse(df$FissionOrBuddingObserved_Genus == 1, 'yes','no')


## ----germline_timing_tidying------------------------------------------------------------------------------------------------------------------------------------
early = c('1','1,2','2')

df$germline_timing_simple<- ifelse(df$germline_timing %in% early, 'early', df$germline_timing) 
df$germline_timing_simple<- ifelse(df$germline_timing == '0', 'no_germline', df$germline_timing_simple) 
df$germline_timing_simple<- ifelse(df$germline_timing == '3', 'adult', df$germline_timing_simple) 


## ----FissionCellNumber------------------------------------------------------------------------------------------------------------------------------------------
#fit the model 
fit_fission_cell_num<-
  brm(data = df,
      family=gaussian(), # family of model
      formula = log(cell_number) ~ 0 + FissionOrBuddingObserved_Genus_nominal,  #formula, the 0 means that there are estimates for both clonal and non-clonal, rather than relative to each other
      iter = 6000000, warmup = 100000, chains = 5,thin = 1000, cores = 5, #chain settings
      prior = prior(normal(0, 10), "b"), #defining the priors- for things in 'class b' set this prior. (can also use get_priors() and set_priors() fucncitons) 
      file = 'fits/fit_FissionCellNumber' ) 

#### Assessing model:  
plot(fit_fission_cell_num) #check that chains converged
pp_check(fit_fission_cell_num) #check the predictions
summary(fit_fission_cell_num) #summary of model 
posterior_summary(fit_fission_cell_num, robust = T)
plot(conditional_effects(fit_fission_cell_num, points = TRUE, ask = F)) 

hyp = hypothesis(fit_fission_cell_num,  "FissionOrBuddingObserved_Genus_nominalno > FissionOrBuddingObserved_Genus_nominalyes") #test hypothesis that there is no difference based on coefficients
hyp
plot(hyp)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------
fit_fission_cell_num$Rhat
autocor(fit_fission_cell_num)


## ----FissionCellType--------------------------------------------------------------------------------------------------------------------------------------------
#fit the model 
fit_fission_cell_type<-
  brm(data = df,
      family=poisson(),
      formula = cell_types ~ 0 + FissionOrBuddingObserved_Genus_nominal + scale(log(cell_number)),  #formula, the 0 means that there are estimates for both clonal and non-clonal, rather than relative to each other
      iter = 6000000, warmup = 100000, chains = 5,thin = 1000, cores = 5, #chain settings
      prior = prior(normal(0, 10), "b"), #defining the priors- for things in 'class b' set this prior. (can also use get_priors() and set_priors() fucncitons) 
      file = 'fits/fit_FissionCellType' ) 

plot(fit_fission_cell_type) #check that chains converged
pp_check(fit_fission_cell_type) #check the predictions
summary(fit_fission_cell_type) #summary of model 
posterior_summary(fit_fission_cell_type, robust = T)
plot(conditional_effects(fit_fission_cell_type, points = TRUE, ask = F)) 

hyp = hypothesis(fit_fission_cell_type,  "FissionOrBuddingObserved_Genus_nominalno > FissionOrBuddingObserved_Genus_nominalyes") #test hypothesis that there is no difference based on coefficients
hyp
plot(hyp)


## ----GermCellNumber---------------------------------------------------------------------------------------------------------------------------------------------
prior <- get_prior(log(cell_number) ~ 0 + germline_timing_simple, family=gaussian(), data = df) #what priors do we need to define?

#fit the model 
fit_germline_cell_num<-
  brm(data = df,
      family=gaussian(), # family of model
      formula = log(cell_number) ~ 0 + germline_timing_simple,  #formula, the 0 means that there are estimates for both clonal and non-clonal, rather than relative to each other
      iter = 6000000, warmup = 100000, chains = 5,thin = 1000, cores = 5, #chain settings
      prior = prior(normal(0, 10), "b"), #defining the priors- for things in 'class b' set this prior. (can also use get_priors() and set_priors() fucncitons) 
      file = 'fits/fit_GermCellNum' ) 

plot(fit_germline_cell_num) #check that chains converged
pp_check(fit_germline_cell_num) #check the predictions
plot(conditional_effects(fit_germline_cell_num), points = TRUE) 

summary(fit_germline_cell_num)
pairs(fit_germline_cell_num)
posterior_summary(fit_germline_cell_num, robust = T) #what are the medians for coefficients? if F, then returns means


## ----GermCellType-----------------------------------------------------------------------------------------------------------------------------------------------
fit_type<-
  brm(data = df,
      family=poisson(),
      formula = cell_types ~ 0 + germline_timing_simple + scale(log(cell_number)),
      iter = 1000000, warmup = 100000, chains = 5, thin = 10000, cores = 5,
      prior = prior(normal(0, 10), "b"), file = 'fits/fit_GermCellType')

plot(fit_type) #check that chains converged
pp_check(fit_type) #check the predictions

summary(fit_type) #summary of model 
posterior_summary(fit_type, robust = T)
plot(conditional_effects(fit_type), points = TRUE, ask = F) 


## ----FissionCellNumber_phy--------------------------------------------------------------------------------------------------------------------------------------
get_prior(data = df,
      family=poisson(), # family of model
      formula = cell_number ~ FissionOrBuddingObserved_Genus_nominal + (1|gr(species_id, cov = CovarMatrix)),
      data2 = list(CovarMatrix = CovarMatrix))

fit_cellnumber_fission_phy<-
  brm(data = df,
      family=poisson(), # family of model
      formula = cell_number ~ FissionOrBuddingObserved_Genus_nominal + (1|gr(species_id, cov = CovarMatrix)),
      prior = c(
        prior(normal(0, 5), "b"),
        prior(normal(0, 20), "Intercept"),
        prior(student_t(3, 0, 20), "sd")), 
      iter = 1000000, warmup = 50000,thin = 100, threads = threading(15), cores = 2, chains = 4,
      file_refit = "always", 
      file = 'fits/fit_phy_FissionCellNumber_long',
      backend = "cmdstanr", 
      data2 = list(CovarMatrix = CovarMatrix))

plot(fit_cellnumber_fission_phy) #check that chains converged
pairs(fit_cellnumber_fission_phy)
pp_check(fit_cellnumber_fission_phy) #check the predictions
summary(fit_cellnumber_fission_phy) #summary of model 
posterior_summary(fit_cellnumber_fission_phy, robust = T)
plot(conditional_effects(fit_cellnumber_fission_phy, points = TRUE, ask = F)) 

hyp = hypothesis(fit_cellnumber_fission_phy,  "FissionOrBuddingObserved_Genus_nominalno > FissionOrBuddingObserved_Genus_nominalyes") #test hypothesis that there is no difference based on coefficients
hyp
plot(hyp)


## ----trying with mcmcglmm---------------------------------------------------------------------------------------------------------------------------------------
# ResolvedPolytomiesTree$node.label<- NULL
# inv.phylo<-inverseA(ResolvedPolytomiesTree,scale=TRUE)
# prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
# 
# df<- df[,-14]
# 
# model_simple<-
#   MCMCglmm(cell_number~FissionOrBuddingObserved_Genus_nominal,
#            random=~ResolvedPolytomiesTree,
#            family="poisson",
#            prior=prior,
#            data=df,
#            nitt=5000,
#            burnin=1000,
#            thin=500)



## ----FissionCellType_phy----------------------------------------------------------------------------------------------------------------------------------------
# #fit the model 
# phylo_fit_fission_cell_type<-
#   brm(data = df,
#       family=poisson(),
#       formula = cell_types ~ 0 + FissionOrBuddingObserved_Genus_nominal + scale(log(cell_number)) + (1|gr(species_id, cov = CovarMatrix)),  #formula, the 0 means that there are estimates for both clonal and non-clonal, rather than relative to each other
#       prior = c(
#         prior(normal(0, 10), "b"),
#         prior(normal(0, 50), "Intercept"),
#         prior(student_t(3, 0, 20), "sd")),  
#       iter = 100000, warmup = 10000, chains = 2, thin = 100, cores = 2, threads = threading(2), 
#       file_refit = "on_change",
#       backend = "cmdstanr", 
#       data2 = list(CovarMatrix = CovarMatrix),
#       file = 'fits/fit_phy_FissionCellType') 
# 
# plot(phylo_fit_fission_cell_type) #check that chains converged
# pp_check(phylo_fit_fission_cell_type) #check the predictions
# summary(phylo_fit_fission_cell_type) #summary of model 
# posterior_summary(phylo_fit_fission_cell_type, robust = T)
# plot(conditional_effects(phylo_fit_fission_cell_type, points = TRUE, ask = F)) 
# 
# hyp = hypothesis(phylo_fit_fission_cell_type,  "FissionOrBuddingObserved_Genus_nominalno > FissionOrBuddingObserved_Genus_nominalyes") #test hypothesis that there is no difference based on coefficients
# hyp
# plot(hyp)


## ----GermCellNum_phy--------------------------------------------------------------------------------------------------------------------------------------------
#fit the model 
# phylofit_germline_cell_num<-
#   brm(data = df,
#       family=poisson(), # family of model
#       formula = cell_number ~ 0 + germline_timing_simple  + (1|gr(species_id, cov = CovarMatrix)),  #formula, the 0 means that there are estimates for both clonal and non-clonal, rather than relative to each other
#       data2 = list(CovarMatrix = CovarMatrix) ,
#       prior = c(
#         prior(normal(0, 10), "b"),
#         prior(normal(0, 50), "Intercept"),
#         prior(student_t(3, 0, 20), "sd")),  
#       iter = 100000, warmup = 10000, chains = 3, thin = 100, cores = 3, threads = threading(6), 
#       backend = "cmdstanr",
#       file_refit = "on_change",
#       file = 'fits/fit_phy_GermCellNum') 
# 
# plot(phylofit_germline_cell_num) #check that chains converged
# pp_check(phylofit_germline_cell_num) #check the predictions
# plot(conditional_effects(phylofit_germline_cell_num), points = TRUE) 
# 
# summary(phylofit_germline_cell_num)
# pairs(phylofit_germline_cell_num)
# posterior_summary(phylofit_germline_cell_num, robust = T) #what are the medians for coefficients? if F, then returns means
# 
# hyp = hypothesis(phylofit_germline_cell_num,  c("germline_timing_simpleearly = germline_timing_simple", "germline_timing_simpleearly = germline_timing_simpleadult", "germline_timing_simpleearly = germline_timing_simpleno_germline")) #test hypothesis that there is no difference based on coefficients
# hyp
# plot(hyp)


## ----GermCellType_phy-------------------------------------------------------------------------------------------------------------------------------------------
# fit_type_phy<-
#   brm(data = df,
#       family=poisson(),
#       formula = cell_types ~ 0 + germline_timing_simple + scale(log(cell_number)) + (1|gr(species_id, cov = CovarMatrix)),
#             prior = c(
#         prior(normal(0, 10), "b"),
#         prior(normal(0, 50), "Intercept"),
#         prior(student_t(3, 0, 20), "sd")),  iter = 100000, warmup = 10000, chains = 3, thin = 100, cores = 3, threads = threading(6), 
#       file_refit = "on_change",
#       backend = "cmdstanr", 
#       file = 'fits/fit_phy_GermCellType', # same simple prior
#       data2 = list(CovarMatrix = CovarMatrix))
# 
# plot(fit_type_phy) #check that chains converged
# pp_check(fit_type_phy) #check the predictions
# summary(fit_type_phy) #summary of model 
# posterior_summary(fit_type_phy, robust = T)
# plot(conditional_effects(fit_type_phy, points = TRUE, ask = F)) 
# 
# hyp = hypothesis(fit_type_phy,  c("germline_timing_simpleearly = germline_timing_simple", "germline_timing_simpleearly = germline_timing_simpleadult", "germline_timing_simpleearly = germline_timing_simpleno_germline"))
# hyp
# plot(hyp)

