#!/usr/bin/env Rscript --vanilla --no-save
library(knitr)
library(argparse)

parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option

parser$add_argument("-i", "--input-file",
                    type = "character",
                    help = "Report Rmd file to knit",
                    metavar = "Rmd report")
parser$add_argument("-o", "--output-file",
                    type = "character",
                    help = "Output filename",
                    default = NULL,
                    metavar = "Rmd report")
  

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

rmarkdown::render(
  input = args$input_file,
  output_format = "pdf_document",
  params = list(skip_LD_api = FALSE),
  output_file = args$output_file
  )
