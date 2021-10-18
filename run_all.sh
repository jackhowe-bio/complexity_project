#!/bin/bash
mkdir -p RScripts/model_outputs
mkdir -p figures
RScript RScripts/OptimisingModels.R
RScript RScripts/ParameterDefinitions.R
RScript RScripts/Model1.R
RScript RScripts/Model2.R
RScript RScripts/Model3.R
RScript RScripts/Model4.R
RScript RScripts/Model5.R
RScript RScripts/Model6.R
RScript RScripts/Model7.R
RScript RScripts/Model8.R

R -e "rmarkdown::render('WorkingNotes.Rmd', output_format='all')"