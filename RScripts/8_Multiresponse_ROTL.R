
df_binary<- df %>%
  mutate(FissionBinary= as.factor(ifelse(df$Fission == 1, 1, 0)), EarlyGermlineBinary= as.factor(ifelse(df$GermNumeric == 1, 1, 0)))



#Here both responses need to 0 and 1s or yes and nos


# run MCMCglmm
  Model_Correlation<- mclapply(1:n_chains, function(i){
   MCMCglmm(cbind(FissionBinary,EarlyGermlineBinary) ~ trait-1,
            random = ~us(trait):species_rotl, #2x2 phylogenetic covariance matrix
            rcov = ~us(trait):units, #2x2 residual covariance matrix
             ginverse=list(species_rotl=inv_tree),family = c("categorical","categorical"),
             data = df_binary,prior=p4, nitt=iterations, burnin=burnin, thin=thinning,verbose = F)
  }, mc.cores = 3)
names(Model_Correlation)<- c('chain1','chain2','chain3','chain4', 'chain5','chain6')

#Model_Correlation<- readRDS(paste(PathForAnalyses, 'R_Objects/ModelOutputs/p4/Model_Correlation_ROTL_p4.RDS', sep = ''))

Model_Correlation_ROTL_Sol<- mcmc.list(lapply(Model_Correlation, function(m) m$Sol))
plot(Model_Correlation_ROTL_Sol) 
gelman.diag(Model_Correlation_ROTL_Sol,multivariate = FALSE)



Model_Correlation_ROTL_VCV<- mcmc.list(lapply(Model_Correlation, function(m) m$VCV))
plot(Model_Correlation_ROTL_VCV) 
gelman.diag(Model_Correlation_ROTL_VCV,multivariate = FALSE)
chain1<- Model_Correlation$chain1
summary(Model_Correlation$chain1)

#*************************************************************************
#CKC edits: 
#phylogenetic heritability in each trait [intraclass correlation coefficient (ICC) for binary traits]
ICC_fission1<-chain1$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']/((chain1$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']+chain1$VCV[,'traitFissionBinary.2:traitFissionBinary.2.units'])+pi^2/3)
posterior.mode(ICC_fission1)
HPDinterval(ICC_fission1)
plot(ICC_fission1)

ICC_germline1<-chain1$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl']/((chain1$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl']+chain1$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.units'])+pi^2/3)
posterior.mode(ICC_germline1)
HPDinterval(ICC_germline1)
plot(ICC_germline1)


#calculate phylogenetic correlation
fission_germline1<-chain1$VCV[,'traitEarlyGermlineBinary.2:traitFissionBinary.2.species_rotl']/sqrt(chain1$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']*chain1$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl'])
plot(fission_germline1)
posterior.mode(fission_germline1)
HPDinterval(fission_germline1)

#Can check the convergence of these terms
chain2<- Model_Correlation$chain2
chain3<- Model_Correlation$chain3
chain4<- Model_Correlation$chain4
chain5<- Model_Correlation$chain5
chain6<- Model_Correlation$chain6

ICC_fission2<-chain2$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']/((chain2$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']+chain2$VCV[,'traitFissionBinary.2:traitFissionBinary.2.units'])+pi^2/3)
ICC_germline2<-chain2$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl']/((chain2$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl']+chain2$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.units'])+pi^2/3)
fission_germline2<-chain2$VCV[,'traitEarlyGermlineBinary.2:traitFissionBinary.2.species_rotl']/sqrt(chain2$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']*chain2$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl'])

ICC_fission3<-chain3$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']/((chain3$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']+chain3$VCV[,'traitFissionBinary.2:traitFissionBinary.2.units'])+pi^2/3)
ICC_germline3<-chain3$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl']/((chain3$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl']+chain3$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.units'])+pi^2/3)
fission_germline3<-chain3$VCV[,'traitEarlyGermlineBinary.2:traitFissionBinary.2.species_rotl']/sqrt(chain3$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']*chain3$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl'])

ICC_fission4<-chain4$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']/((chain4$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']+chain4$VCV[,'traitFissionBinary.2:traitFissionBinary.2.units'])+pi^2/3)
ICC_germline4<-chain4$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl']/((chain4$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl']+chain4$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.units'])+pi^2/3)
fission_germline4<-chain4$VCV[,'traitEarlyGermlineBinary.2:traitFissionBinary.2.species_rotl']/sqrt(chain4$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']*chain4$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl'])

ICC_fission5<-chain5$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']/((chain5$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']+chain5$VCV[,'traitFissionBinary.2:traitFissionBinary.2.units'])+pi^2/3)
ICC_germline5<-chain5$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl']/((chain5$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl']+chain5$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.units'])+pi^2/3)
fission_germline5<-chain5$VCV[,'traitEarlyGermlineBinary.2:traitFissionBinary.2.species_rotl']/sqrt(chain5$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']*chain5$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl'])

ICC_fission6<-chain6$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']/((chain6$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']+chain6$VCV[,'traitFissionBinary.2:traitFissionBinary.2.units'])+pi^2/3)
ICC_germline6<-chain6$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl']/((chain6$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl']+chain6$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.units'])+pi^2/3)
fission_germline6<-chain6$VCV[,'traitEarlyGermlineBinary.2:traitFissionBinary.2.species_rotl']/sqrt(chain6$VCV[,'traitFissionBinary.2:traitFissionBinary.2.species_rotl']*chain6$VCV[,'traitEarlyGermlineBinary.2:traitEarlyGermlineBinary.2.species_rotl'])

ICC_fission<-mcmc.list(ICC_fission1,ICC_fission2,ICC_fission3,ICC_fission4,ICC_fission5,ICC_fission6)
ICC_germline<-mcmc.list(ICC_germline1,ICC_germline2,ICC_germline3,ICC_germline4,ICC_germline5,ICC_germline6)
fission_germline<-mcmc.list(fission_germline1,fission_germline2,fission_germline3,fission_germline4,fission_germline5,fission_germline6)
varnames(fission_germline) <- 'correlation'


plot(ICC_fission)
gelman.diag(ICC_fission,multivariate = FALSE)

plot(ICC_germline)
gelman.diag(ICC_germline,multivariate = FALSE)

plot(fission_germline)
gelman.diag(fission_germline,multivariate = FALSE)


HPDinterval(fission_germline)


summary(fission_germline)

trace<- mcmc_trace(fission_germline) + theme_minimal()
areas<- mcmc_areas(fission_germline) + theme_minimal()
intervals<- mcmc_intervals(fission_germline) + theme_minimal()

a<- ggarrange(trace, areas, intervals, ncol = 3, labels = 'AUTO')
a


pdf(paste(PathForAnalyses, 'R_Objects/ModelOutputs/p4/Fission_Germline.pdf', sep = ''), width = 30, height = 10)
print(a)
dev.off()

#*************************************************************************
 Model_Correlation_Sol_gg<- ggmcmc::ggs(Model_Correlation_ROTL_Sol)
 
 diagnostic_filename<- paste(PathForAnalyses, "R_Objects/ModelDiagnostics/p4/Model_Correlation_ROTL_p4.pdf", sep = '')
 ggmcmc::ggmcmc(Model_Correlation_Sol_gg, file = diagnostic_filename)

MCMCglmm_filename<- paste(PathForAnalyses, 'R_Objects/ModelOutputs/p4/Model_Correlation_ROTL_p4.RDS', sep = '')
saveRDS(Model_Correlation, MCMCglmm_filename)
 
 

# 
# #Univariate regression is covariance between response and predictor/variance in predictor
# Model_Correlation_VCV<- mcmc.list(lapply(Model_Correlation, function(m) m$VCV))
# plot(Model_Correlation_VCV) 
# gelman.diag(Model_Correlation_VCV,multivariate = FALSE)
# 
# Model_Correlation_VCV_gg<- ggmcmc::ggs(Model_Correlation_VCV)
# 
# # for between species use the parameters ending in species 
# Covariance<- Model_Correlation_VCV_gg %>%
#   filter(Parameter == 'traitFissionBinary.2:traitEarlyGermlineBinary.2.species_rotl') %>% 
#   rename(Covariance = value) %>% select(-Parameter)
# 
# 
# #in previous model was: 'traitFissionBinary.1:traitEarlyGermlineBinary.1.species'
# 
# VariancePredictor<- Model_Correlation_VCV_gg %>%
#   filter(Parameter == 'traitFissionBinary.2:traitFissionBinary.2.species_rotl') %>%
#   rename(PredictorVariance = value) %>% select(-Parameter)
# 
# # in previous model was: traitFissionBinary.1:traitFissionBinary.1.species
# 
# Regression<- merge(Covariance, VariancePredictor) %>%
#   mutate(value = Covariance / PredictorVariance)
# glimpse(Regression)
# 
# pdf('RScripts/ModelOutputs/p4/Regression.pdf')
# ggplot(Regression, aes(x = Iteration, y = value, colour = as.factor(Chain))) + geom_line(size = 1, alpha = 0.5) + geom_smooth() + theme_minimal() + ylab('Regression coefficient')
# ggplot(Regression, aes(x = Iteration, y = Covariance, colour = as.factor(Chain))) + geom_line(size = 1, alpha = 0.5) + geom_smooth() + theme_minimal()
# ggplot(Regression, aes(x = Iteration, y = PredictorVariance, colour = as.factor(Chain))) + geom_line(size = 1, alpha = 0.5) + geom_smooth() + theme_minimal()
# 
# dev.off()
# 
# ggplot(Regression, aes(x = Iteration, y = value, colour = as.factor(Chain))) + geom_line(size = 1, alpha = 0.5) + geom_smooth() + theme_minimal() + ylab('Regression coefficient')
