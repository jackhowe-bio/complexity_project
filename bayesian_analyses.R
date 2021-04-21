# attempting analysis

#load packages
library(ape)
library(ggplot2)
library(brms) #https://rdrr.io/cran/brms/f/vignettes/brms_phylogenetics.Rmd
library(rotl) #see https://cran.r-project.org/web/packages/rotl/vignettes/rotl.html


df<- read.csv('germline_data_1.0.csv')

# trying to recreate Berti's data: does clonality lead to more cells?


prior <- get_prior(cell_number ~ clonal, family = bernoulli(link = "logit"), data = df)
prior$prior[3]


fit_clonal<-
  brm(data = df,
      family = gaussian(),
      formula = cell_number ~ clonal, 
      iter = 6000000, warmup = 1000000, chains = 3, seed = 8, thin = 1000,
      prior = prior,
      file = "fits/fit1")

plot(fit_clonal)
summary(fit_clonal)

pairs(fit_clonal,
      off_diag_args = list(size = 1/3, alpha = 1/3))

hyp = hypothesis(fit_clonal, "clonalclonal > clonalnonMclonal")

plot(hyp)
