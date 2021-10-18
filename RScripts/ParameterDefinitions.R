iterations<- 8000000
burnin<- 1000000
thinning <- 100
n_chains<- 6

saveRDS(iterations, file = 'RScripts/iterations.RDS')
saveRDS(burnin, file = 'RScripts/burnin.RDS')
saveRDS(thinning, file = 'RScripts/thinning.RDS')
saveRDS(n_chains, file = 'RScripts/n_chains.RDS')
