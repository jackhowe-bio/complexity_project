# attempting analysis

#load packages
library(ape)
library(ggplot2)
library(tidyverse)
library(brms) #https://rdrr.io/cran/brms/f/vignettes/brms_phylogenetics.Rmd
library(tidybayes)
library(rotl) #see https://cran.r-project.org/web/packages/rotl/vignettes/rotl.html

setwd('/Users/user/Documents/oxford/comparative_data/tidy/analyses/complexity_project')
df<- read.csv('germline_data_1.0.csv')

# trying to recreate Berti's data: does clonality lead to more obligate?
df %>%
  group_by(clonal, Obligate_or_facultative) %>%
  count() #summary table just to check that the results kinda fit

prior <- get_prior(Obligate_or_facultative ~ 0 + clonal, family=bernoulli(link = 'logit'), data = df) #defining what priors we need
prior$prior[2] <- "normal(0, 1)" #default priors are flat, changed them here to a normal distribution
prior$prior[3] <- "normal(0, 1)"

#fit the model 
fit_obl<-
  brm(data = df,
      family=bernoulli(link = 'logit'),
      formula = Obligate_or_facultative ~ 0 + clonal, 
      iter = 6000000, warmup = 1000000, chains = 5, seed = 8, thin = 1000, cores = 5,
      prior = prior)


plot(fit_obl) #check that chains converged

pp_check(fit_obl) #check the predictions

summary(fit_obl) #summary of model 
posterior_summary(fit_obl, robust = T) #what are the medians for coefficients?

post <- posterior_samples(fit_obl) #take the posterior samples
post <- #convert them back to probabilities (from log odds)
  post %>% 
  summarise(theta_clonal = exp(b_clonalclonal) / (1 + exp(b_clonalclonal)),
         theta_nonclonal     = exp(b_clonalnonMclonal)/ (1 + exp(b_clonalnonMclonal))) %>% 
  mutate(`theta_clonal - theta_nonclonal` = theta_clonal - theta_nonclonal)

gathered_post <- # pivot table to make plot 
  post %>% 
  select(starts_with("theta")) %>% 
  gather() %>% 
  mutate(key = factor(key, levels = c("theta_clonal", "theta_nonclonal", "theta_clonal - theta_nonclonal")))

gathered_post %>% #draw plot
  ggplot(aes(x = value, y = 0, fill = key)) +
  stat_histinterval(point_interval = mode_hdi, .width = .95,
                    slab_color = "white", outline_bars = T,
                    normalize = "panels") +
  scale_y_continuous(NULL, breaks = NULL) +
  theme_minimal() +
  facet_wrap(~key, scales = 'free') 

hyp = hypothesis(fit_obl,  "clonalclonal = clonalnonMclonal") #test hypothesis that there is no difference based on coefficients
hyp
plot(hyp)



