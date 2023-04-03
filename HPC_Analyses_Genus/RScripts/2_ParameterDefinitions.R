
#setwd
setwd('/Users/pcx971/Documents/oxford/complexity/complexity_project/HPC_Analyses/')

# define priors to be used
p1=list(R = list(V = 1, nu=0.002), 
        G = list(G1=list(V = 1, nu = 0.002)))

saveRDS(p1, file = 'RScripts/R_Objects/priors_1.RDS')

p2=list(R=list(V=1, nu=1), 
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.var=1000)))
saveRDS(p2, file = 'RScripts/R_Objects/priors_2.RDS')

p3=list(R=list(V=1, nu=2),
        G=list(G1=list(V=1, nu=2, alpha.mu=0, alpha.var=1000)))
saveRDS(p3, file = 'RScripts/R_Objects/priors_3.RDS')

#setting iterations etc.
iterations<- 8000000
burnin<- 1000000
thinning <- 1000
n_chains<- 6


# for testing
#iterations<- 10000
#burnin<- 1000
#thinning <- 10
#n_chains <- 6

saveRDS(iterations, file = 'RScripts/R_Objects/iterations.RDS')
saveRDS(burnin, file = 'RScripts/R_Objects/burnin.RDS')
saveRDS(thinning, file = 'RScripts/R_Objects/thinning.RDS')
saveRDS(n_chains, file = 'RScripts/R_Objects/n_chains.RDS')
