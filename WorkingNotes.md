*Working Notes:* Single-celled bottlenecks, germlines and the evolution
of complex multi-cellularity
================





-   [Data Collection](#data-collection)
    -   [Search strategy](#search-strategy)
    -   [Production of phylogenetic
        tree](#production-of-phylogenetic-tree)
-   [Statistical Analyses](#statistical-analyses)
    -   [MCMCglmm parameters](#mcmcglmm-parameters)
    -   [Without Phylogeny](#without-phylogeny)
        -   [**Model 1**: Fission vs Cell
            Number](#model-1-fission-vs-cell-number)
        -   [**Model 2**: Fission vs Cell
            Types](#model-2-fission-vs-cell-types)
        -   [**Model 3**: Germline vs Cell
            Numbers](#model-3-germline-vs-cell-numbers)
        -   [**Model 4**: Germline vs Cell
            Types](#model-4-germline-vs-cell-types)
    -   [Phylogenetically Informed
        Models](#phylogenetically-informed-models)
        -   [**Model 5**: Fission vs Cell
            Number](#model-5-fission-vs-cell-number)
        -   [**Model 6**: Fission vs Cell
            Types](#model-6-fission-vs-cell-types)
        -   [**Model 7**: Germline vs Cell
            Number](#model-7-germline-vs-cell-number)
        -   [**Model 8**: Germline vs Cell
            Types](#model-8-germline-vs-cell-types)
    -   [Open Questions](#open-questions)
-   [References](#references)

# Data Collection

## Search strategy

Searches were conducted broadly for literature focussed on reproductive
mode and germline development across the tree of life. This included
chapters reviews and chapters within textbooks.

As Fisher paper was used for estimates of individual complexity (which
in turn used Bell paper as a foundation), narrower searches were
conducted for each species/genus from Bell paper on Web of Knowledge and
on Google Scholar.

`(ALL = (reproduct* OR sex* OR asex* OR vegetat* OR fissi* OR clonal* OR regenerat* OR rhizo* OR germ-line* OR germline* OR germ line* OR bud* OR fragment* OR parthenogen* OR stolon*)) AND (ALL = TAXON)`

We recorded whether sexual, parthenogenetic/clonal and agametic have
been observed as binary values. We did not attempt to capture the
relative frequency of different reproductive strategies, as these data
do not exist for the majority of species.

Research into reproduction and development is heterogeneous across
organisms, and this is necessarily reflected in the required searching
effort for different groups: the reproductive biology and development of
model organisms such as C. elegans, M. musculus, D. melanogaster, A.
thaliana are well known, but in many other groups reproduction may never
have been observed. We therefore conducted searches for each species,
but if no literature discussing reproductive strategy was observed, then
we conducted an additional search at the genus level. This assumes that
genera will tend to be relatively similar in their reproductive
strategies: this is not always the case. The planarian S. mediterranea,
for example, has strictly sexual and strictly asexual strains within
even the same species. However, we are focussed on patterns through
longer spans of evolution than between individual genera.

Caveats:

-   Sexual organisms contain more cells because they have gonads that
    asexual organisms lack…  
-   hard to fit algae with gametophyte/sporophyte stages: they have life
    cycles where some stages can fragment, where some stages can
    reproduce by spores, parthenogenesis or by sex. If an algae can stay
    in one loop and reproduce exclusively through
    parthenogenesis/fragmentation, then counts– similar to organisms
    that can fission, but don’t always.

## Production of phylogenetic tree

We used a phylogenetic tree to control for non-independence of species
based on shared evolutionary history. The tree was constructed within R
using the latest evolutionary classifications found on the Tree of Life,
AlgaeBase.org, and the World Register of Marine species. The
relationships among species were reconstructed by ordering the taxa from
Kingdom through to species, and grouping according to these names.

As a comparison, we also constructed a tree using the ‘R Tree of Life
Project.’ These two trees were largely congruent: some larger groups had
switched places, but within these groups relationships were
predominantly the same. As the Rtol tree dropped X data points from the
tree, we used the tree based on the taxa names. Multichotomies within
the tree were randomly resolved, before branch lengths were generated as
described by \[@grafen1990\]. Branches smaller than 10^{-25} were
deleted, and the dichotomies here collapsed to multichotomies. Figure ()
shows a cophylogeny based on each tree.

# Statistical Analyses

## MCMCglmm parameters

All analyses were conducted in R \[@R-base\] using the package MCMCglmm
\[@MCMCglmm\], while documents were produced using \[RMarkdown\]. All
data and code are accessible at github.

Model parameters were optimised using the first model described in our
results, for which we ran a total of 38 MCMCglmm chains of varying
lengths (500000 - 10000000 iterations), with varying warm-ups (100000 -
1000000, and with thinning of either 100 or 1000 fold, see Figure
S@ref(fig:OptimisationFigure). All subsequent models were then fit using
the combination of these parameters where the autocorrelation of
successive sampled mean and variance were minimal: 8^{6} iterations, a
warm-up of 10^{6} iterations and thinning by a factor of 100. In all
fitted models, the autocorrelation was well below the suggested
tolerable maximum of 0.1 \[@hadfield?\]. For each model, 6 chains were
run which were visually inspected for chain convergence. Convergence was
also supported by the Gelman-Rubin \[@Gelman-Rubin\] convergence
diagnostic, which approximated 1 (\<1.05) in all cases–these are
reported in the summary of each model below.

The four models described in section @ref(Without Phylogeny) are
phylogenetically naive, and treat each species as independent data
points. The default priors used for fixed effects, and residual variance
prior of V = 1 and nu = 0.002.

The differences between each level for the fixed effects were calculated
at each MCMC iteration to produce a posterior distribution for the
difference. Levels are considered statistically significant if the 95%
credible interval of this difference distribution did not overlap with
0, and if the proportion of MCMC iterations that were greater or less
than 0 was less than 0.05.

The entire analysis, including the parameter optimisation step and
creating all output documents, runs in approximately 3hrs 15 mins using
a 2020 MacBook Pro running 4 chains in parallel.

![Autocorrelation of successively sampled mean and variance values from
posterior
distribution](WorkingNotes_files/figure-gfm/OptimisationFigure-1.png)

## Without Phylogeny

### **Model 1**: Fission vs Cell Number

![**Model 1: Cell Numbers vs Fission** *A* Traceplots for the estimated
means, *B* Estimates for means from posterior distribution, dots
represent median, thick and thin lines indicate 90% and 95% of highest
posterior density regions, respectively. *C* Density plot of estimated
differences between fissioning and non-fissioning organisms, bar
represents 90% and 95% credible
intervals.](WorkingNotes_files/figure-gfm/Model1Mean-1.png)

| Fixed Effects                     | Posterior Mode (CI) | pMCMC |
|:----------------------------------|:--------------------|------:|
| FissionOrBuddingObserved_Species0 | 7.02 (5.18, 9.03)   | 0.000 |
| FissionOrBuddingObserved_Species1 | 13.01 (9.59, 15.58) | 0.000 |
| FissionOrBuddingObserved_Species? | 8.84 (-3.36, 21.12) | 0.140 |
| FissionOrBuddingObserved_Species  | 6.91 (1.23, 12.28)  | 0.018 |

Table S1.1: Estimates of Fixed Effects

| Fixed Effects Comparisons                                              | Posterior Mode (CI)   | pMCMC |
|:-----------------------------------------------------------------------|:----------------------|------:|
| FissionOrBuddingObserved_Species0 vs FissionOrBuddingObserved_Species1 | -5.62 (-8.97, -1.83)  | 0.003 |
| FissionOrBuddingObserved_Species0 vs FissionOrBuddingObserved_Species? | -1.31 (-14.33, 10.47) | 0.738 |
| FissionOrBuddingObserved_Species0 vs FissionOrBuddingObserved_Species  | 1.37 (-5.6, 6.12)     | 0.894 |
| FissionOrBuddingObserved_Species1 vs FissionOrBuddingObserved_Species? | 4.44 (-9.21, 15.99)   | 0.592 |
| FissionOrBuddingObserved_Species1 vs FissionOrBuddingObserved_Species  | 5.96 (-0.48, 12.05)   | 0.067 |
| FissionOrBuddingObserved_Species? vs FissionOrBuddingObserved_Species  | 2.84 (-11.23, 15.83)  | 0.718 |

Table S1.2: Comparison Between Fixed Effects

**Do organisms that reproduce by fission have more cells?** Fissiparous
organisms appear to be larger (makes sense, trees, fungi, algae, etc)

*priors*: p1=list(R = list(V = 1, nu=0.002)) #sets prior for residual
variance, the defaults are used as priors for fixed effects (see
MCMCglmm course notes)

### **Model 2**: Fission vs Cell Types

*priors* p1=list(R = list(V = 1, nu=0.002))

![**Model 2: Cell Types vs Fission** *A* Traceplots for the estimated
means, *B* Estimates for means from posterior distribution, dots
represent median, thick and thin lines indicate 90% and 95% of highest
posterior density regions, respectively. *C* Density plot of estimated
differences between fissioning and non-fissioning organisms, bar
represents 90% and 95% credible
intervals.](WorkingNotes_files/figure-gfm/Model2Mean-1.png)

| Fixed Effects                     | Posterior Mode (CI) | pMCMC |
|:----------------------------------|:--------------------|------:|
| FissionOrBuddingObserved_Species0 | 1.86 (1.65, 2.05)   | 0.000 |
| FissionOrBuddingObserved_Species1 | 1.98 (1.7, 2.21)    | 0.000 |
| FissionOrBuddingObserved_Species? | 2.51 (1.49, 3.44)   | 0.000 |
| FissionOrBuddingObserved_Species  | 1.5 (0.73, 2.4)     | 0.001 |
| scale(log(cell_number))           | 0.8 (0.65, 0.97)    | 0.000 |

Table S1.1: Estimates of Fixed Effects

| Fixed Effects Comparisons                                              | Posterior Mode (CI) | pMCMC |
|:-----------------------------------------------------------------------|:--------------------|------:|
| FissionOrBuddingObserved_Species0 vs FissionOrBuddingObserved_Species1 | -0.09 (-0.42, 0.22) | 0.532 |
| FissionOrBuddingObserved_Species0 vs FissionOrBuddingObserved_Species? | -0.73 (-1.59, 0.41) | 0.236 |
| FissionOrBuddingObserved_Species0 vs FissionOrBuddingObserved_Species  | 0.32 (-0.58, 1.15)  | 0.517 |
| FissionOrBuddingObserved_Species0 vs scale(log(cell_number))           | 1.09 (0.79, 1.32)   | 0.000 |
| FissionOrBuddingObserved_Species1 vs FissionOrBuddingObserved_Species? | -0.58 (-1.5, 0.5)   | 0.331 |
| FissionOrBuddingObserved_Species1 vs FissionOrBuddingObserved_Species  | 0.33 (-0.5, 1.26)   | 0.390 |
| FissionOrBuddingObserved_Species1 vs scale(log(cell_number))           | 1.21 (0.82, 1.47)   | 0.000 |
| FissionOrBuddingObserved_Species? vs FissionOrBuddingObserved_Species  | 0.86 (-0.41, 2.17)  | 0.178 |
| FissionOrBuddingObserved_Species? vs scale(log(cell_number))           | 1.7 (0.64, 2.66)    | 0.001 |
| FissionOrBuddingObserved_Species vs scale(log(cell_number))            | 0.72 (-0.07, 1.6)   | 0.077 |

Table S1.2: Comparison Between Fixed Effects

**Do organisms that reproduce by fission have more cell types?** HCI
overlaps with zero, so doesn’t seem likely.

### **Model 3**: Germline vs Cell Numbers

Should we subset to only those organisms that have sterile cells for the
germline models?

![**Model 3: Cell number vs Germline** *A* Traceplots for the estimated
means, *B* Estimates for means from posterior distribution, dots
represent median, thick and thin lines indicate 90% and 95% of highest
posterior density regions, respectively. *C* Density plot of estimated
differences between early germline segregators and organisms that
segregate germlines continuously as adults, bar represents 90% and 95%
credible intervals.](WorkingNotes_files/figure-gfm/Model3Mean-1.png)

| Fixed Effects                     | Posterior Mode (CI)  | pMCMC |
|:----------------------------------|:---------------------|------:|
| germline_timing_simple3           | 5.38 (-7.18, 20.88)  | 0.334 |
| germline_timing_simpleadult       | 16.71 (15.4, 18.17)  | 0.000 |
| germline_timing_simpleearly       | 14.73 (11.67, 17.34) | 0.000 |
| germline_timing_simpleno_germline | 3.15 (0.26, 6.11)    | 0.033 |

Table S1.1: Estimates of Fixed Effects

| Fixed Effects Comparisons                                        | Posterior Mode (CI)  | pMCMC |
|:-----------------------------------------------------------------|:---------------------|------:|
| germline_timing_simple3 vs germline_timing_simpleadult           | -9.61 (-24.21, 4.01) | 0.170 |
| germline_timing_simple3 vs germline_timing_simpleearly           | -6.61 (-21.82, 6.75) | 0.295 |
| germline_timing_simple3 vs germline_timing_simpleno_germline     | 3.99 (-10.69, 17.92) | 0.610 |
| germline_timing_simpleadult vs germline_timing_simpleearly       | 2.22 (-0.91, 5.39)   | 0.166 |
| germline_timing_simpleadult vs germline_timing_simpleno_germline | 13.74 (10.31, 16.76) | 0.000 |
| germline_timing_simpleearly vs germline_timing_simpleno_germline | 11.33 (7.33, 15.48)  | 0.000 |

Table S1.2: Comparison Between Fixed Effects

**Do organisms with early segregating germline have more cells?** HCI
just about overlaps with 0, so maayyyybe, but not clear.

*priors* p1=list(R = list(V = 1, nu=0.002))

### **Model 4**: Germline vs Cell Types

![**Model 4: Cell Types vs Germline** *A* Traceplots for the estimated
means, *B* Estimates for means from posterior distribution, dots
represent median, thick and thin lines indicate 90% and 95% of highest
posterior density regions, respectively. *C* Density plot of estimated
differences between early germline segregators and organisms that
segregate germlines continuously as adults, bar represents 90% and 95%
credible intervals.](WorkingNotes_files/figure-gfm/Model4Mean-1.png)

| Fixed Effects                     | Posterior Mode (CI) | pMCMC |
|:----------------------------------|:--------------------|------:|
| germline_timing_simpleadult       | 1 (0.68, 1.3)       | 0.000 |
| germline_timing_simpleearly       | 1.53 (0.87, 2.13)   | 0.000 |
| germline_timing_simpleno_germline | 0.59 (0.12, 1.11)   | 0.019 |
| scale(log(cell_number))           | 0.73 (0.52, 0.98)   | 0.000 |

Table S1.1: Estimates of Fixed Effects

| Fixed Effects Comparisons                                        | Posterior Mode (CI) | pMCMC |
|:-----------------------------------------------------------------|:--------------------|------:|
| germline_timing_simpleadult vs germline_timing_simpleearly       | -0.58 (-1.23, 0.16) | 0.136 |
| germline_timing_simpleadult vs germline_timing_simpleno_germline | 0.38 (-0.28, 1.02)  | 0.278 |
| germline_timing_simpleadult vs scale(log(cell_number))           | 0.23 (-0.27, 0.72)  | 0.354 |
| germline_timing_simpleearly vs germline_timing_simpleno_germline | 0.85 (0.1, 1.69)    | 0.031 |
| germline_timing_simpleearly vs scale(log(cell_number))           | 0.85 (0.08, 1.42)   | 0.032 |
| germline_timing_simpleno_germline vs scale(log(cell_number))     | -0.09 (-0.59, 0.33) | 0.603 |

Table S1.2: Comparison Between Fixed Effects

**Do organisms that segregate germline early have more cell types?**
Again, just about overlaps with 0, so not clear.

p1=list(R = list(V = 1, nu=0.002))

## Phylogenetically Informed Models

Models below here use inverse covariance matrix describing the
relationships among species to control for phylogeny.

### **Model 5**: Fission vs Cell Number

![**Model 5: Cell Number vs Fission with phylogeny ** *A* Traceplots for
the estimated means, *B* Estimates for means from posterior
distribution, dots represent median, thick and thin lines indicate 90%
and 95% of highest posterior density regions, respectively. *C* Density
plot of estimated differences between fissioning and non-fissioning
organisms, bar represents 90% and 95% credible
intervals.](WorkingNotes_files/figure-gfm/Model5Mean-1.png)

Again, just about overlaps with 0, so not clear.

| Fixed Effects                     | Posterior Mode (CI) | pMCMC |
|:----------------------------------|:--------------------|------:|
| FissionOrBuddingObserved_Species0 | 8.92 (0.38, 18.32)  | 0.044 |
| FissionOrBuddingObserved_Species1 | 9.99 (1.73, 19.62)  | 0.020 |
| FissionOrBuddingObserved_Species? | 8.44 (-0.63, 18.63) | 0.061 |
| FissionOrBuddingObserved_Species  | 9.85 (0.69, 19.19)  | 0.040 |

Table S1.1: Estimates of Fixed Effects

| Fixed Effects Comparisons                                              | Posterior Mode (CI) | pMCMC |
|:-----------------------------------------------------------------------|:--------------------|------:|
| FissionOrBuddingObserved_Species0 vs FissionOrBuddingObserved_Species1 | -1.54 (-3.98, 1)    | 0.223 |
| FissionOrBuddingObserved_Species0 vs FissionOrBuddingObserved_Species? | 0.01 (-3.76, 3.86)  | 0.999 |
| FissionOrBuddingObserved_Species0 vs FissionOrBuddingObserved_Species  | -0.43 (-3.63, 2.63) | 0.751 |
| FissionOrBuddingObserved_Species1 vs FissionOrBuddingObserved_Species? | 1.37 (-3.1, 5.96)   | 0.500 |
| FissionOrBuddingObserved_Species1 vs FissionOrBuddingObserved_Species  | 0.92 (-1.81, 3.88)  | 0.465 |
| FissionOrBuddingObserved_Species? vs FissionOrBuddingObserved_Species  | -0.35 (-5.28, 4.55) | 0.839 |

Table S1.2: Comparison Between Fixed Effects

p2=list(R = list(V = 1, nu=0.002), G = list(G1=list(V=1, nu=0.002)))

### **Model 6**: Fission vs Cell Types

![**Model 6: Cell Number vs Fission with phylogeny ** *A* Traceplots for
the estimated means, *B* Estimates for means from posterior
distribution, dots represent median, thick and thin lines indicate 90%
and 95% of highest posterior density regions, respectively. *C* Density
plot of estimated differences between fissioning and non-fissioning
organisms, bar represents 90% and 95% credible
intervals.](WorkingNotes_files/figure-gfm/Model6Mean-1.png)

Overlaps with 0, no difference

p2=list(R = list(V = 1, nu=0.002), G = list(G1=list(V=1, nu=0.002)))

| Fixed Effects                     | Posterior Mode (CI) | pMCMC |
|:----------------------------------|:--------------------|------:|
| FissionOrBuddingObserved_Species0 | 1.07 (-0.01, 2.09)  | 0.049 |
| FissionOrBuddingObserved_Species1 | 1.14 (0.17, 2.26)   | 0.029 |
| FissionOrBuddingObserved_Species? | 1.05 (-0.42, 2.36)  | 0.164 |
| FissionOrBuddingObserved_Species  | 0.84 (-0.43, 1.99)  | 0.189 |
| scale(log(cell_number))           | 0.65 (0.34, 0.93)   | 0.000 |

Table S1.1: Estimates of Fixed Effects

| Fixed Effects Comparisons                                              | Posterior Mode (CI) | pMCMC |
|:-----------------------------------------------------------------------|:--------------------|------:|
| FissionOrBuddingObserved_Species0 vs FissionOrBuddingObserved_Species1 | -0.17 (-0.63, 0.37) | 0.599 |
| FissionOrBuddingObserved_Species0 vs FissionOrBuddingObserved_Species? | 0.05 (-0.85, 1.08)  | 0.815 |
| FissionOrBuddingObserved_Species0 vs FissionOrBuddingObserved_Species  | 0.33 (-0.5, 1.04)   | 0.493 |
| FissionOrBuddingObserved_Species0 vs scale(log(cell_number))           | 0.38 (-0.68, 1.54)  | 0.418 |
| FissionOrBuddingObserved_Species1 vs FissionOrBuddingObserved_Species? | 0.3 (-0.82, 1.29)   | 0.654 |
| FissionOrBuddingObserved_Species1 vs FissionOrBuddingObserved_Species  | 0.39 (-0.31, 1.1)   | 0.259 |
| FissionOrBuddingObserved_Species1 vs scale(log(cell_number))           | 0.64 (-0.53, 1.71)  | 0.300 |
| FissionOrBuddingObserved_Species? vs FissionOrBuddingObserved_Species  | 0.18 (-1.05, 1.35)  | 0.771 |
| FissionOrBuddingObserved_Species? vs scale(log(cell_number))           | 0.39 (-1.07, 1.79)  | 0.612 |
| FissionOrBuddingObserved_Species vs scale(log(cell_number))            | 0.1 (-1.13, 1.44)   | 0.781 |

Table S1.2: Comparison Between Fixed Effects

### **Model 7**: Germline vs Cell Number

![**Model 7: Cell Number vs Fission with phylogeny ** *A* Traceplots for
the estimated means, *B* Estimates for means from posterior
distribution, dots represent median, thick and thin lines indicate 90%
and 95% of highest posterior density regions, respectively. *C* Density
plot of estimated differences between early germline segregators and
organisms that segregate germlines continuously as adults, bar
represents 90% and 95% credible
intervals.](WorkingNotes_files/figure-gfm/Model7Mean-1.png)

| Fixed Effects                     | Posterior Mode (CI)  | pMCMC |
|:----------------------------------|:---------------------|------:|
| germline_timing_simple3           | 8.49 (-15.17, 28.24) | 0.545 |
| germline_timing_simpleadult       | 14.32 (-3.62, 34.17) | 0.106 |
| germline_timing_simpleearly       | 11.94 (-8.16, 30.07) | 0.226 |
| germline_timing_simpleno_germline | 12.1 (-6.96, 31.31)  | 0.204 |

Table S1.1: Estimates of Fixed Effects

| Fixed Effects Comparisons                                        | Posterior Mode (CI)  | pMCMC |
|:-----------------------------------------------------------------|:---------------------|------:|
| germline_timing_simple3 vs germline_timing_simpleadult           | -8.7 (-19.97, 2.02)  | 0.113 |
| germline_timing_simple3 vs germline_timing_simpleearly           | -5.63 (-16.1, 6.14)  | 0.378 |
| germline_timing_simple3 vs germline_timing_simpleno_germline     | -5.67 (-17.09, 6.06) | 0.343 |
| germline_timing_simpleadult vs germline_timing_simpleearly       | 3.76 (0.42, 7.48)    | 0.031 |
| germline_timing_simpleadult vs germline_timing_simpleno_germline | 3.71 (0.01, 6.72)    | 0.048 |
| germline_timing_simpleearly vs germline_timing_simpleno_germline | -0.22 (-5.43, 4.4)   | 0.823 |

Table S1.2: Comparison Between Fixed Effects

Things with a germline might be smaller: the 95% CI is *just* below 0
(-0.34)

p2=list(R = list(V = 1, nu=0.002), G = list(G1=list(V=1, nu=0.002)))

### **Model 8**: Germline vs Cell Types

    ## 
    ##  Iterations = 1000001:7999901
    ##  Thinning interval  = 100
    ##  Sample size  = 70000 
    ## 
    ##  DIC: 226.7579 
    ## 
    ##  G-structure:  ~species
    ## 
    ##         post.mean l-95% CI u-95% CI eff.samp
    ## species    0.7402   0.2043    1.418    34743
    ## 
    ##  R-structure:  ~units
    ## 
    ##       post.mean  l-95% CI u-95% CI eff.samp
    ## units    0.0123 0.0001529  0.04578    50977
    ## 
    ##  Location effects: cell_types ~ germline_timing_simple - 1 + scale(log(cell_number)) 
    ## 
    ##                                   post.mean l-95% CI u-95% CI eff.samp    pMCMC
    ## germline_timing_simpleadult          1.1211   0.2479   2.0623    55252  0.02006
    ## germline_timing_simpleearly          1.7797   0.7297   2.8027    50297  0.00271
    ## germline_timing_simpleno_germline    0.7392  -0.3472   1.7870    30346  0.16497
    ## scale(log(cell_number))              0.5993   0.3087   0.8953    35619 8.57e-05
    ##                                      
    ## germline_timing_simpleadult       *  
    ## germline_timing_simpleearly       ** 
    ## germline_timing_simpleno_germline    
    ## scale(log(cell_number))           ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

![**Model 8: Cell Number vs Fission with phylogeny ** *A* Traceplots for
the estimated means, *B* Estimates for means from posterior
distribution, dots represent median, thick and thin lines indicate 90%
and 95% of highest posterior density regions, respectively. *C* Density
plot of estimated differences between early germline segregators and
organisms that segregate germlines continuously as adults, bar
represents 90% and 95% credible
intervals.](WorkingNotes_files/figure-gfm/Model8Mean-1.png)

Seems like they may be smaller, but that they have more cell types per
cell– HCI is just above 0 (0.0379847).

p2=list(R = list(V = 1, nu=0.002), G = list(G1=list(V=1, nu=0.002)))

Is pMCMC just the number of simulated cases where difference is \<0? In
which case:

| Fixed Effects                     | Posterior Mode (CI) | pMCMC |
|:----------------------------------|:--------------------|------:|
| germline_timing_simpleadult       | 1.13 (0.25, 2.06)   | 0.020 |
| germline_timing_simpleearly       | 1.78 (0.73, 2.8)    | 0.003 |
| germline_timing_simpleno_germline | 0.64 (-0.35, 1.79)  | 0.165 |
| scale(log(cell_number))           | 0.57 (0.31, 0.9)    | 0.000 |

Table S1.1: Estimates of Fixed Effects

| Fixed Effects Comparisons                                        | Posterior Mode (CI)  | pMCMC |
|:-----------------------------------------------------------------|:---------------------|------:|
| germline_timing_simpleadult vs germline_timing_simpleearly       | -0.69 (-1.26, -0.04) | 0.037 |
| germline_timing_simpleadult vs germline_timing_simpleno_germline | 0.37 (-0.35, 1.1)    | 0.301 |
| germline_timing_simpleadult vs scale(log(cell_number))           | 0.54 (-0.47, 1.54)   | 0.292 |
| germline_timing_simpleearly vs germline_timing_simpleno_germline | 1 (0.15, 1.96)       | 0.026 |
| germline_timing_simpleearly vs scale(log(cell_number))           | 1.22 (0.05, 2.25)    | 0.042 |
| germline_timing_simpleno_germline vs scale(log(cell_number))     | 0.09 (-0.93, 1.23)   | 0.786 |

Table S1.2: Comparison Between Fixed Effects

## Open Questions

Phylogenetic correlation between germline and fission– multivariate
model? How to test this?

# References
