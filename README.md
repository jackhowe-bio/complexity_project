# Single-celled bottlenecks, germlines and the evolution of complex multi-cellularity

This repository contains all data and code for analysis.

## Running the files

Cloning this github repository should allow the script (germline_bayesian_analyses.Rmd) to be run relatively quickly (if packages installed). 

The fisher_bayesian_analyses.rmd file replicates the results from Berti's paper (with shorter chains), as a test that things generally work as they should (they do!). The germline_bayesian_analyses.rmd runs the analyses using the new dataset. They currently use the genus-level info for reproduction, the analyses on the species-level data are next on the list!

Output as germline_bayesian_analyses.md, which can be read pretty well online here. 

The input data is saved in the 'data' directory, and the chains from the brms runs are saved in the 'fits' directory. The 'brm' functions that run the analyses contain a "file = 'fits/...'" parameter, which tells it to load the run that is saved in the 'fits' folder and ignore any other parameters. Removing the fit parameter, or changing the file name will cause the analysis to be rerun (or using an 'update.brmsfit' function to update the fit might be quicker?)

## Directory structure:
```
-- README.md (this file)  
-- fisher_bayesian_analyses.Rmd (Rmd that produces results from Fisher paper) 
-- germline_bayesian_analyses.Rmd (the working R script for the novel analyses) 
-- germline_bayesian_analyses.md (a readable output of the script with figures, etc. ) 
-- data/ (sub dir)  
  |-- germline_data_1.1.csv (csv containing dataset)  
  |-- phylogeny_all.txt (phylogeny produced by Rmd script)  
  |-- phylogeny_all_res_polytomy.txt (same phylogeny, but with polytomies resolved randomly)  
  |-- phylogeny_species_names.csv (produced to match species names in dataset to the labels used by rotl phylogeny)
  |-- phylogeny_species_names_in_tree.csv (same dictionary, but subset to include only those that are present in the rotl tree)  
-- fits/ (sub dir that contains the sampling data from MCMC chains, means that the models do not need re-run every time) BUT, if parameters changed, then this file needs deleted, or changed within code)
-- germline_bayesian_analyses_files/ (directory that contains the images used in the .md output, ignore this folder)
```

## Search strategy

Searches were conducted broadly for literature focussed on reproductive mode and germline development across the tree of life. This included chapters reviews and chapters within textbooks.

As Fisher paper was used for estimates of individual complexity (which in turn used Bell paper as a foundation), narrower searches were conducted for each species/genus from Bell paper on Web of Knowledge and on Google Scholar. 

`(ALL = (reproduct* OR sex* OR asex* OR vegetat* OR fissi* OR clonal* OR regenerat* OR rhizo* OR germ-line* OR germline* OR germ line* OR bud* OR fragment* OR parthenogen* OR stolon*)) AND (ALL = TAXON)` 

We recorded whether sexual, parthenogenetic/clonal and agametic have been observed as binary values. We did not attempt to capture the relative frequency of different reproductive strategies, as these data do not exist for the majority of species. 

Research into reproduction and development is heterogeneous across organisms, and this is necessarily reflected in the required searching effort for different groups: the reproductive biology and development of model organisms such as C. elegans, M. musculus, D. melanogaster, A. thaliana is well known, but in many other groups reproduction may never have been observed. We therefore conducted searches for each species, but if no literature discussing reproductive strategy was observed, then we conducted an additional search at the genus level. This assumes that genera will tend to be relatively similar in their reproductive strategies: this is not always the case. The planarian S. mediterranea, for example, has strictly sexual and strictly asexual strains within even the same species. However, we are focussed on patterns through longer spans of evolution than between individual genera. 

Caveats:

* Sexual organisms contain more cells because they have gonads that asexual organisms lack...  
* hard to fit algae with gametophyte/sporophyte stages: they have life cycles where some stages can fragment, where some stages can reproduce by spores, parthenogenesis or by sex. If an algae can stay in one loop and reproduce exclusively through parthenogenesis/fragmentation, then counts-- similar to organisms that can fission, but don't always.


