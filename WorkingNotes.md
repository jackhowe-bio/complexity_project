*Working Notes:* Single-celled bottlenecks, germlines and the evolution
of complex multi-cellularity
================
truetruetrue

# Data Collection

Data were collected by…

Phylogeny was constructed by

# Statistical Analyses

All analyses were conducted in R (R Core Team 2021) using the package
MCMCglmm \[MCMCglmm\], while documents were produced using
\[RMarkdown\]. All data and code are accessible at \[github\].

Model parameters were optimised using the first model described in our
results, for which we ran a total of 38 MCMCglmm chains of varying
lengths (500000 - 10000000 iterations), with varying warm-ups (100000 -
1000000, and with thinning of either 100 or 1000 fold, see Figure
S@ref(fig:OptimisationFigure). All subsequent models were then fit using
the combination of these parameters where the autocorrelation in the
sampled mean and variance were minimal: \[X\] iterations, a warm-up of
\[Y\] iterations and thinning by a factor of \[Z\]. In all fitted
models, the autocorrelation was well below the suggested tolerable
maximum of 0.1 (**hadfield?**). For each model, \[n_chains\] chains were
run. These were visually inspected for chain convergence, and the
Gelman-Rubin (**Gelman-Rubin?**) convergence diagnostic approximated 1
(\<1.05) in all cases–these are reported in the summary of each model
below.

![Fig S1: Autocorrelation of successively sampled mean and variance
values from posterior
distribution](WorkingNotes_files/figure-gfm/OptimisationFigure-1.png) ##
Testing for a correlation between reproductive biology and multicellular
complexity

## Phylogenetically Informed Tests

# References

<div id="refs" class="references csl-bib-body hanging-indent"
custom-style="Bibliography">

<div id="ref-R-base" class="csl-entry">

R Core Team. 2021. *R: A Language and Environment for Statistical
Computing*. Vienna, Austria: R Foundation for Statistical Computing.
<https://www.R-project.org/>.

</div>

</div>
