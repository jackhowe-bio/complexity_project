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
Project’. These two trees were largely congruent: some larger groups had
switched places, but within these groups relationships were
predominantly the same. As the Rtol tree dropped X data points from the
tree, we used the tree based on the taxa names. Multichotomies within
the tree were randomly resolved, before branch lengths were generated as
described by (**grafen1989?**). Branches smaller than 10^{-25} were
deleted, and the dichotomies here collapsed to multichotomies. Figure ()
shows a cophylogeny based on each tree.

# Statistical Analyses

## MCMCglmm parameters

All analyses were conducted in R (R Core Team 2021) using the package
MCMCglmm (**MCMCglmm?**), while documents were produced using
\[RMarkdown\]. All data and code are accessible at github..

Model parameters were optimised using the first model described in our
results, for which we ran a total of 38 MCMCglmm chains of varying
lengths (500000 - 10000000 iterations), with varying warm-ups (100000 -
1000000, and with thinning of either 100 or 1000 fold, see Figure
S@ref(fig:OptimisationFigure). All subsequent models were then fit using
the combination of these parameters where the autocorrelation of
successive sampled mean and variance were minimal: 8^{6} iterations, a
warm-up of 10^{6} iterations and thinning by a factor of 100. In all
fitted models, the autocorrelation was well below the suggested
tolerable maximum of 0.1 (**hadfield?**). For each model, 6 chains were
run which were visually inspected for chain convergence. Convergence was
also supported by the Gelman-Rubin (**Gelman-Rubin?**) convergence
diagnostic, which approximated 1 (\<1.05) in all cases–these are
reported in the summary of each model below.

![Autocorrelation of successively sampled mean and variance values from
posterior
distribution](WorkingNotes_files/figure-gfm/OptimisationFigure-1.png)

## Without Phylogeny

### **Model 1**: Fission vs Cell Number

![**Model 1: Cell Numbers vs Fission** *A* Traceplots for the estimated
means f, *B* Estimates for means from posterior distribution, dots
represent median, thick and thin lines indicate 90% and 95% of highest
posterior density regions,
respectively.](WorkingNotes_files/figure-gfm/Model1Mean-1.png)

### **Model 2**: Fission vs Cell Types

![**Model 2: Cell Types vs Fission** *A* Traceplots for the estimated
means f, *B* Estimates for means from posterior distribution, dots
represent median, thick and thin lines indicate 90% and 95% of highest
posterior density regions,
respectively.](WorkingNotes_files/figure-gfm/Model2Mean-1.png)

### **Model 3**: Germline vs Cell Numbers

Should we subset to only those organisms that have sterile cells for the
germline models?

![**Model 3: Cell number vs Germline** *A* Traceplots for the estimated
means f, *B* Estimates for means from posterior distribution, dots
represent median, thick and thin lines indicate 90% and 95% of highest
posterior density regions,
respectively.](WorkingNotes_files/figure-gfm/Model3Mean-1.png)

### **Model 4**: Germline vs Cell Types

![**Model 4: Cell Types vs Germline** *A* Traceplots for the estimated
means f, *B* Estimates for means from posterior distribution, dots
represent median, thick and thin lines indicate 90% and 95% of highest
posterior density regions,
respectively.](WorkingNotes_files/figure-gfm/Model4Mean-1.png)

## Phylogenetically Informed Models

The datapoints are not independent: they have shared evolutionary
history of varying degrees. Should we exclude some of the

### **Model 5**: Fission vs Cell Number

![**Model 5: Cell Number vs Fission with germline ** *A* Traceplots for
the estimated means f, *B* Estimates for means from posterior
distribution, dots represent median, thick and thin lines indicate 90%
and 95% of highest posterior density regions,
respectively.](WorkingNotes_files/figure-gfm/Model5Mean-1.png)

### **Model 6**: Fission vs Cell Types

![**Model 6: Cell Number vs Fission with germline ** *A* Traceplots for
the estimated means f, *B* Estimates for means from posterior
distribution, dots represent median, thick and thin lines indicate 90%
and 95% of highest posterior density regions,
respectively.](WorkingNotes_files/figure-gfm/Model6Mean-1.png)

### **Model 7**: Germline vs Cell Number

![**Model 7: Cell Number vs Fission with germline ** *A* Traceplots for
the estimated means f, *B* Estimates for means from posterior
distribution, dots represent median, thick and thin lines indicate 90%
and 95% of highest posterior density regions,
respectively.](WorkingNotes_files/figure-gfm/Model7Mean-1.png)

### **Model 8**: Germline vs Cell Types

![**Model 8: Cell Number vs Fission with germline ** *A* Traceplots for
the estimated means f, *B* Estimates for means from posterior
distribution, dots represent median, thick and thin lines indicate 90%
and 95% of highest posterior density regions,
respectively.](WorkingNotes_files/figure-gfm/Model8Mean-1.png) ## Open
Questions

Phylogenetic correlation between germline and fission– multivariate
model? How to test this?

# References

<div id="refs" class="references csl-bib-body hanging-indent"
custom-style="Bibliography">

<div id="ref-R-base" class="csl-entry">

R Core Team. 2021. *R: A Language and Environment for Statistical
Computing*. Vienna, Austria: R Foundation for Statistical Computing.
<https://www.R-project.org/>.

</div>

</div>
