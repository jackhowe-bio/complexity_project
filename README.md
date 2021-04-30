# Single-celled bottlenecks, germlines and the evolution of complex multi-cellularity

This repository contains all data and code for analysis.

## Directory structure:
```
-- README.md (this file)  
-- fisher_bayesian_analyses.Rmd (Rmd that produces results from Fisher paper) 
-- germline_bayesian_analyses.Rmd (the working R script for the novel analyses) 
-- data/ (sub dir)  
  |-- germline_data_1.0.csv (csv containing dataset)  
  |-- phylogeny_all.txt (phylogeny produced by Rmd script)  
  |-- phylogeny_all_res_polytomy.txt (same phylogeny, but with polytomies resolved randomly)  
  |-- phylogeny_species_names.csv (produced to match species names in dataset to the labels used by rotl phylogeny)
  |-- phylogeny_species_names_in_tree.csv (same dictionary, but subset to include only those that are present in the rotl tree)  
-- fits/ (sub dir that contains the sampling data from MCMC chains, means that the models do not need re-run every time) BUT, if parameters changed, then this file needs deleted, or changed within code)  
```
