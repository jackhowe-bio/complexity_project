#!/bin/bash

## commented out the MCMCglmm scripts once they've been run to save from running them again

RScript --verbose 1_DataPrep.R # not working because of rotl server?
RScript --verbose 2_ParameterDefinitions.R
RScript --verbose 3_Optimisation.R
RScript --verbose 4_Model1_ROTL.R
RScript --verbose 5_Model2_ROTL.R
RScript --verbose 6_Model3_ROTL.R
RScript --verbose 7_Model4_ROTL.R
RScript --verbose 8_Model1And2_Interpretation.R
RScript --verbose 9_Model3And4_Interpretation.R
RScript --verbose 10_Figures.R
RScript --verbose 11_Multiresponse_ROTL.R
RScript --verbose 12_GermlineWithoutCellNumber.R