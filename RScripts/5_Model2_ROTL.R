# Model 5: fission vs number of cells / phylogeny
#Read in the priors


i = 1
for(prior_set in list(p1,p2,p3)){
Model2_ROTL_parallel <- mclapply(1:n_chains, function(i){
  MCMCglmm(Types ~ Fission-1 + scale(log(Number)), #-1 here removes the intercept equivalent to 0 in brms
           random = ~species_rotl, ginverse=list(species_rotl=inv_tree), # phylogeny modelled by linking species to inverse distance matrix created from phylogeny
           family ="poisson",data = df,prior=prior_set, nitt=iterations, burnin=burnin, thin=thinning ,verbose = F, pr=F)
}, mc.cores = 3)
names(Model2_ROTL_parallel)<- c('chain1','chain2','chain3','chain4', 'chain5','chain6')


Model2_ROTL_Sol<- mcmc.list(lapply(Model2_ROTL_parallel, function(m) m$Sol))
#plot(Model2_ROTL_Sol) 
gelman.diag(Model2_ROTL_Sol,multivariate = FALSE)

Model2_ROTL_VCV<- mcmc.list(lapply(Model2_ROTL_parallel, function(m) m$VCV))
#plot(Model2_ROTL_VCV) 
gelman.diag(Model2_ROTL_VCV,multivariate = FALSE)
chain1<- Model2_ROTL_parallel$chain1
summary(Model2_ROTL_parallel$chain1)

M2_Sol_gg<- ggmcmc::ggs(Model2_ROTL_Sol)

diagnostic_filename<- paste(PathForAnalyses, 'R_Objects/ModelDiagnostics/p', i,'/Model2_ROTL_p',i, '.pdf', sep = '')
ggmcmc::ggmcmc(M2_Sol_gg, file = diagnostic_filename)


MCMCglmm_filename<- paste(PathForAnalyses, 'R_Objects/ModelOutputs/p', i,'/Model2_ROTL_p',i, '.RDS', sep = '')
saveRDS(Model2_ROTL_parallel, MCMCglmm_filename)
i = i + 1
####################
}