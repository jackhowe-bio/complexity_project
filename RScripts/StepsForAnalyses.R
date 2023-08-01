# then prepare the dataset and the tree: this saves to the file that is loaded next
source("RScripts/1_DataPrep.R", verbose = T)

#Read in the data and the tree
File=paste(PathForAnalyses, 'R_Objects/metadata.RDS', sep = '')
df = readRDS(File)
File=paste(PathForAnalyses, 'R_Objects/inv_tree.RDS', sep = '')
inv_tree = readRDS(File)

# run the models
#source("RScripts/3_Optimisation.R", verbose = T)
source("RScripts/4_Model1_ROTL.R")
source("RScripts/5_Model2_ROTL.R")
source("RScripts/6_Model3_ROTL.R")
source('RScripts/7_Model4_ROTL.R')

# then the correlation models
source('RScripts/8_Multiresponse_ROTL.R')
source('RScripts/9_GermlineWithoutCellNumber.R')
source('RScripts/10_AncestralStateReconstructions.R')

# draw the figures
source('RScripts/11_Figures.R')


