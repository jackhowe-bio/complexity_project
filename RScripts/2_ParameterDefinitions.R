# define priors to be used
p1=list(R = list(V = 1, nu=0.002), 
        G = list(G1=list(V = 1, nu = 0.002)))

File <- paste(PathForAnalyses, 'R_Objects/priors_1.RDS', sep = "")
saveRDS(p1, file = File)

p2=list(R=list(V=1, nu=1), 
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.var=1000)))
File <- paste(PathForAnalyses, 'R_Objects/priors_2.RDS', sep = "")
saveRDS(p2, file = File)

p3=list(R=list(V=1, nu=2),
        G=list(G1=list(V=1, nu=2, alpha.mu=0, alpha.var=1000)))
File <- paste(PathForAnalyses, 'R_Objects/priors_3.RDS', sep = "")
saveRDS(p3, file = File)



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
File <- paste(PathForAnalyses, 'R_Objects/iterations.RDS', sep = "")
saveRDS(iterations, file = 'RScripts/R_Objects/iterations.RDS')

File <- paste(PathForAnalyses, 'R_Objects/burnin.RDS', sep = "")
saveRDS(burnin, file = 'RScripts/R_Objects/burnin.RDS')

File <- paste(PathForAnalyses, 'R_Objects/priors_3.RDS', sep = "")
saveRDS(thinning, file = 'RScripts/R_Objects/thinning.RDS')

File <- paste(PathForAnalyses, 'R_Objects/priors_3.RDS', sep = "")
saveRDS(n_chains, file = 'RScripts/R_Objects/n_chains.RDS')
