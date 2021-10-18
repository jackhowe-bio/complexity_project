# Single-celled bottlenecks, germlines and the evolution of complex multi-cellularity


## Running the files

This repository contains all data and code for analysis.

The analyses is described in full in **WorkingNotes.Rmd** and its output **WorkingNotes.pdf**. The Rmd file contains snippets of code for creating figures and reporting results, but all analyses are conducted within smaller scripts for each model. 


## Directory structure:

```
.
├── README.md
├── RScripts
│   ├── DataImportPhylogenyConstruction.R
│   ├── Model1.R
│   ├── Model2.R
│   ├── Model3.R
│   ├── Model4.R
│   ├── Model5.R
│   ├── Model6.R
│   ├── Model7.R
│   ├── Model8.R
│   ├── OptimisingModels.R
│   ├── OptimisingModels.Rdata  (not uploaded)
│   ├── ParameterDefinitions.R
│   ├── burnin.RDS
│   ├── iterations.RDS
│   ├── model_outputs (not uploaded)
│   │   ├── Model1.RDS (not uploaded)
│   │   ├── Model2.RDS (not uploaded)
│   │   ├── Model3.RDS (not uploaded)
│   │   ├── Model4.RDS (not uploaded)
│   │   ├── Model5.RDS (not uploaded)
│   │   ├── Model6.RDS (not uploaded)
│   │   ├── Model7.RDS (not uploaded)
│   │   └── Model8.RDS (not uploaded)
│   ├── n_chains.RDS
│   └── thinning.RDS
├── WorkingNotes.Rmd
├── WorkingNotes_files
│   └── figure-latex
│       ├── Model1Mean-1.pdf
│       ├── Model2Mean-1.pdf
│       ├── Model3Mean-1.pdf
│       ├── Model4Mean-1.pdf
│       ├── Model5Mean-1.pdf
│       ├── Model6Mean-1.pdf
│       ├── Model7Mean-1.pdf
│       └── OptimisationFigure-1.pdf
├── complexity_project.Rproj
├── data
│   ├── References\ 1.docx
│   ├── germline_data_1.2.csv
│   └── germline_data_1.2.numbers
└── figures
    └── phylogenies.pdf
```

The .pdf and .md reports should be readable without running any scripts, but to create the outputs from the files in order to run the script locally, run the following command: 

```
bash run_all.sh
```

This runs a shell script that creates the folders required for outputs, runs all models and then knits the output from the RMarkdown Document 'WorkingNotes.Rmd' to a html document, a pdf document, and a github friendly document (when closer to submission will output a word document)


