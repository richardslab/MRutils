
library(knitr); 
library(argparse)

parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument("-t", "--token", 
                    type = "string",
                    help = "LD-Link API access token. If you don't have one go here: https://ldlink.nci.nih.gov/?tab=apiaccess",
                    metavar = "token")
parser$add_argument("-i", "--input-file", 
                    type = "string", 
                    help = "Report Rmd file to knit", 
                    metavar = "Rmd report")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

print(args)
rmarkdown::render(args.input_file, 'pdf_document', params = list(LDLink_token = args.token))
