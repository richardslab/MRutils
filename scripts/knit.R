#!env Rscript 

library(knitr); 
setwd("scripts/")
rmarkdown::render('VD_COVID_MR_ALL.Rmd', 'pdf_document')
