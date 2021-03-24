
library(knitr); 
args = commandArgs(trailingOnly=TRUE)
print(args)
rmarkdown::render(args[1], 'pdf_document')
