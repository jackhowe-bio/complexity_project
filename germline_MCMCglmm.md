Single-celled bottlenecks, germlines and the evolution of complex
multi-cellularity
================

-   [Pre-amble](#pre-amble)
    -   [load packages](#load-packages)
    -   [read in data](#read-in-data)
    -   [data tidying](#data-tidying)
    -   [create phylogeny](#create-phylogeny)
-   [Analyses](#analyses)
    -   [Phylogenetically naive](#phylogenetically-naive)

This document follows the same structure as
‘germline_bayesian_analyses.Rmd’, but uses MCMCglmm package rather than
BRMS as brms was slow to run and difficult to get to converge.

# Pre-amble

``` r
knitr::opts_chunk$set(message=F) # set multiple global options here
```

## load packages

``` r
library(ape)
library(tidyverse)
library(MCMCglmm)
```

## read in data

``` r
df<- read_csv('data/germline_data_1.2.csv', col_types = list( #import data with defining column types
  cell_types = col_integer(),
  cell_number = col_integer(),
  clonal = col_factor(),
  SexualReproductionObserved_Species = col_factor(),
  AsexualGameticObserved_Species = col_factor(),
  FissionOrBuddingObserved_Species = col_factor(),
  SexualReproductionObserved_Genus = col_factor(),
  AsexualGameticObserved_Genus = col_factor(),
  FissionOrBuddingObserved_Genus = col_factor(),
  fragmentation_or_budding = col_factor(),
  germline_timing = col_factor()
))

#create a genus column
df<- df %>% 
  separate(species, into = c('genus', 'species_name'), remove = F) 
```

    ## Warning: Expected 2 pieces. Additional pieces discarded in 3 rows [167, 168,
    ## 169].

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 16 rows [3, 10,
    ## 35, 39, 41, 66, 84, 102, 103, 109, 110, 111, 134, 135, 140, 153].

``` r
glimpse(df)
```

    ## Rows: 196
    ## Columns: 28

    ## Warning: One or more parsing issues, see `problems()` for details

    ## $ species                            <chr> "Acrasis rosea", "Schimmelmannia da…
    ## $ genus                              <chr> "Acrasis", "Schimmelmannia", "Acyto…
    ## $ species_name                       <chr> "rosea", "dawsonii", NA, "pallidum"…
    ## $ species_original_name              <chr> "Acrasis rosea", "Schimmelmannia da…
    ## $ Obligate_or_facultative            <chr> "facultative", "obligate", "faculta…
    ## $ Simple_or_complex                  <chr> "simple", "complex", "simple", "sim…
    ## $ Sterile_soma                       <chr> "no", "yes", "no", "yes", "yes", "y…
    ## $ Proportion_of_sterile_cells        <dbl> 0, NA, 0, 45, NA, NA, NA, NA, NA, N…
    ## $ cell_types                         <int> 2, 11, 1, 2, 12, 14, 42, 9, 12, 3, …
    ## $ cell_number                        <int> NA, NA, NA, NA, NA, NA, NA, NA, NA,…
    ## $ clonal                             <fct> non-clonal, clonal, non-clonal, non…
    ## $ `Fisher Ref`                       <chr> "[18-20]", "[21]", "[20, 47]", "[48…
    ## $ kingdom                            <chr> "Protozoa", "Plantae", "Protozoa", …
    ## $ phylum                             <chr> "Percolozoa", "Rhodophyta", "Myceto…
    ## $ class                              <chr> "Heterolobosea", "Florideophyceae",…
    ## $ order                              <chr> "Acrasida", "Acrosymphytales", "Acy…
    ## $ family                             <chr> "Acrasiaceae", "Acrosymphytaceae", …
    ## $ SexualReproductionObserved_Species <fct> 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0,…
    ## $ AsexualGameticObserved_Species     <fct> 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0,…
    ## $ FissionOrBuddingObserved_Species   <fct> 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1,…
    ## $ SexualReproductionObserved_Genus   <fct> 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1,…
    ## $ AsexualGameticObserved_Genus       <fct> 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,…
    ## $ FissionOrBuddingObserved_Genus     <fct> 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1,…
    ## $ notes                              <chr> "Life cycle in Brown & Silberman", …
    ## $ fragmentation_or_budding           <fct> NA, , f, f, fission, , b, NA, b, NA…
    ## $ germline_timing                    <fct> "3", "3", "3", "3", "1,2", "3", "3"…
    ## $ refs                               <chr> "-77", "(2, 30)", "-77", "-78", "(1…
    ## $ ...26                              <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA,…

``` r
df$cell_number<- as.integer(as.numeric(df$cell_number))
```

## data tidying

``` r
early = c('1','1,2','2')

df$germline_timing_simple<- ifelse(df$germline_timing %in% early, 'early', df$germline_timing)  # early germline segregation
df$germline_timing_simple<- ifelse(df$germline_timing == '0', 'no_germline', df$germline_timing_simple) # no germline/soma divide
df$germline_timing_simple<- ifelse(df$germline_timing == '3', 'adult', df$germline_timing_simple) # germline segregated in adult tissues
```

## create phylogeny

``` r
#Get species in dataset and order by classification
species<-df %>% select(kingdom, phylum, class, order,family,genus,species) %>% arrange(kingdom, phylum, class, order,family,genus,species)

species<- species %>% mutate(across(.cols = everything(),~as.factor(.x))) #make factors
tree<-as.phylo(~kingdom/phylum/class/order/family/genus/species, data=species) #draw tree
```

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in f.rec(subtaxo[subtaxo[, u] == levels[l], ][1:(u - 1)]): leaves names
    ## are not unique.

    ## Warning in as.phylo.formula(~kingdom/phylum/class/order/family/genus/species, :
    ## NAs introduced by coercion

``` r
tree<-multi2di(tree, random=T) #collapse and resolve multichotomies (randomly)
tree<-makeNodeLabel(tree, method="number", prefix="Node") #make node labels
tree<-compute.brlen(tree, method = "Grafen", power = 1) #compute branch lengths
tree<-di2multi(tree, tol=1e-25) # collapse and resolve multichotomies (again?)
tree<-chronoMPL(tree) #estimates age of nodes from distance to all tips
treeTrimmed<- drop.tip(tree, which(is.na(tree$tip.label))) #drop any tips where the label is NA (not sure this is an issue any longer)

#should we include an outgroup? 

is.binary.tree(treeTrimmed)
```

    ## [1] TRUE

``` r
is.ultrametric(treeTrimmed)
```

    ## [1] TRUE

``` r
pdf('figures/phylogenies.pdf')
plot(treeTrimmed, cex = 0.4)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

# Analyses

Conducting Bayesian analyses using MCMCglmm, first without phylogeny
included, then including phylogeny.

Each analysis uses relatively uninformative priors, 5 chains, 6,000,000
iterations, 100000 of which are discarded as warm-up, thinned by a
factor of 1000.

(Currently use reproductive traits of genus rather than species)

-   Phylogenetically naive:
    -   Does a single-celled bottleneck correlate with increased cell
        number? (priors = normal dist, mean of 0, sd of 10)
    -   Does a single-celled bottleneck correlate with increased cell
        types (per cell)? prior(normal(0, 10), “b”)
    -   Does germline timing correlate with increased cell number? prior
        = prior(normal(0, 10), “b”)
    -   Does germline timing correlate with increased cell types (per
        cell)?
-   Phylogenetically informed:
    -   Does a single-celled bottleneck correlate with increased cell
        number?  
    -   Does a single-celled bottleneck correlate with increased cell
        types (per cell)?
    -   Does germline timing correlate with increased cell number?
    -   Does germline timing correlate with increased cell types (per
        cell)?

To do: do any of these things need scaled?

## Phylogenetically naive

#CKC edit: MCMCglmm models ##Without phylogeny
