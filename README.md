# Single-celled bottlenecks, germlines and the evolution of complex multi-cellularity


## Running the files

This repo contains all data and scripts for all analyses run.

Most things are now stored in the HPC_Analyses folder (named because I had ambitions 
of running on the server, but had issues with R libraries there so ran locally).

All scripts are in the `RScripts` folder. 



## Directory structure:
```
.
|-- HPC_Analyses
|   |-- Data
|   |   |-- AllMetadata.txt
|   |   |-- Phylogeny_ROTL.txt
|   |   |-- Phylogeny_ROTL_intermediate.txt
|   |   |-- Phylogeny_ROTL_original.txt
|   |   |-- ROTL_metadata.txt
|   |   `-- germline_data_1.4.csv
|   |-- RScripts
|   |   |-- 10_Figures.R #generates some of the figures for the manuscript (except anvio fig)
|   |   |-- 11_Multiresponse_ROTL.R
|   |   |-- 1_DataPrep.R #generates the phylogeny and preps the data
|   |   |-- 2_ParameterDefinitions.R #defines iterations etc, needs edited after running 3
|   |   |-- 3_Optimisation.R #runs combinations of mcmcglmm iterations/burnin/thinning for testing
|   |   |-- 4_Model1_ROTL.R 
|   |   |-- 5_Model2_ROTL.R
|   |   |-- 6_Model3_ROTL.R
|   |   |-- 7_Model4_ROTL.R
|   |   |-- 8_Model1And2_Interpretation.R 
|   |   |-- 9_Model3And4_Interpretation.R
|   |   |-- ModelDiagnostics #contains pdfs with diagnostics from mcmcglmm
|   |   |   |-- p1
|   |   |   |   `-- Model3_ROTL_p1.pdf
|   |   |   |-- p2
|   |   |   |   `-- Model3_ROTL_p2.pdf
|   |   |   |-- p3
|   |   |   |   `-- Model3_ROTL_p3.pdf
|   |   |   `-- p4
|   |   |       `-- Model_Correlation_ROTL_p4.pdf
|   |   |-- ModelOutputs #the RDS files are large so are not uploaded
|   |   |   |-- p1 
|   |   |   |   |-- Model1SummaryFigure_p1.pdf #pdfs show the HPD intervals and chains
|   |   |   |   |-- Model1_ROTL_p1.RDS 
|   |   |   |   |-- Model2SummaryFigure_p1.pdf
|   |   |   |   |-- Model2_ROTL_p1.RDS
|   |   |   |   |-- Model3SummaryFigure_p1.pdf
|   |   |   |   |-- Model3_ROTL_p1.RDS
|   |   |   |   |-- Model4SummaryFigure_p1.pdf
|   |   |   |   `-- Model4_ROTL_p1.RDS
|   |   |   |-- p2
|   |   |   |   |-- Model1SummaryFigure_p2.pdf
|   |   |   |   |-- Model1_ROTL_p2.RDS
|   |   |   |   |-- Model2SummaryFigure_p2.pdf
|   |   |   |   |-- Model2_ROTL_p2.RDS
|   |   |   |   |-- Model3SummaryFigure_p2.pdf
|   |   |   |   |-- Model3_ROTL_p2.RDS
|   |   |   |   |-- Model4SummaryFigure_p2.pdf
|   |   |   |   `-- Model4_ROTL_p2.RDS
|   |   |   |-- p3
|   |   |   |   |-- Model1SummaryFigure_p3.pdf
|   |   |   |   |-- Model1_ROTL_p3.RDS
|   |   |   |   |-- Model2SummaryFigure_p3.pdf
|   |   |   |   |-- Model2_ROTL_p3.RDS
|   |   |   |   |-- Model3SummaryFigure_p3.pdf
|   |   |   |   |-- Model3_ROTL_p3.RDS
|   |   |   |   |-- Model4SummaryFigure_p3.pdf
|   |   |   |   `-- Model4_ROTL_p3.RDS
|   |   |   `-- p4
|   |   |       |-- Model_Correlation_ROTL_p4.RDS
|   |   |       `-- Regression.pdf
|   |   |-- R_Objects #contains priors and parameters for mcmcglmm
|   |   |   |-- burnin.RDS
|   |   |   |-- inv_tree.RDS
|   |   |   |-- iterations.RDS
|   |   |   |-- metadata.RDS
|   |   |   |-- n_chains.RDS
|   |   |   |-- priors_1.RDS
|   |   |   |-- priors_2.RDS
|   |   |   |-- priors_3.RDS
|   |   |   `-- thinning.RDS
|   |   `-- run_all.sh # runs each of the scripts above (needs updated)
|   |-- anvio
|   |   |-- AnvioOutput.svg #output from anvio
|   |   |-- FigureProfile.db #needed for anvio file
|   |   |-- Phylogeny_ROTL.txt #generated from data prep
|   |   `-- metadata.txt #generated from dataprep
|   `-- figures #figures for the manuscript
|       |-- Figure2.pdf
|       |-- FigureBars.pdf
|       |-- FigureFission.pdf
|       |-- FigureGerm.pdf
|       |-- Optimisation.pdf
|       `-- phylogeniesROTL.pdf
|-- MCMCglmm\ functions.R #not sure if needed
|-- README.html 
|-- README.md
|-- complexity_project.Rproj
|-- data #used previously, tidier data in HPC folder
`-- manuscript


```

