# Conflict Reducing Innovations in Development Enable Increased Multicellular Complexity


## README

This repo contains all data and scripts for all analyses run for the manuscript. 

The `data` folder contains the data and associated references. 

There are then 2 directories: 

- `HPC_Analyses` contains the data with genus level data excluded
- `HPC_Analyses_Genus` contains the data with genus level data included

Within each of these directories are the data, scripts and outputs.

- `anvio` contains the inputs and outputs for the anvio software used to draw and label 
the phylogeny
- `Data` includes the data files
- `figures` contains the output summary figures used in the manuscript
- `RScripts` contains all the scripts used to generate results and figures
  - .R files starting with 1 through 12 run all the analyses. To run all files, 
  use the `run_all.sh` script while in this folder. 
  - `SupplementaryInfo.Rmd` and `SupplementaryInfo.pdf` are in the `HPC_Analyses` 
  folder. This generates the supplementary info. (a copy is present in the first dir)
  -`ModelOutputs` contains the figures that were used to populate the supplementary 
  info file for each prior set. The .RDS model objects generated in the R Scripts
  are not uploaded here because of their large size. They are generated when the scripts are run. 


## Running the analyses

When the dir is downloaded, run the following to conduct all analyses and produce all figures:

```{bash}
cd HPC_Analyses/RScripts
bash run_all.sh

cd ..
cd ..
cd HPC_Analyses_Genus/RScripts
bash run_all.sh
```

This may involve the installation of the relavant R packages if they are not already installed.

Chains for each model are run in parallel to speed up analyses (3 at a time). This can be sped up if your machine allows by changing the `mc.cores = 3` line to a number >3. 

Once run the `SupplementaryInfo.Rmd` file can be knit to create the supplement. 
