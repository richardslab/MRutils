#!/usr/bin/env Rscript

### select high significant snps from study:

#libraries
library(tidyr)
library(readxl)

#parameters
p_threshold <- 5e-8
comman_threshold <- 0.01
skip_rows <- 4
input <- "1-s2.0-S0002929720300173-mmc2.xlsx"
output <- "significant_VD_SNPs.tsv"
columns <- c(1:5, 13:17)

# read the data with no name repair (since the columns have repreating names)
# subset to the study of choice
data <- read_excel(input, skip = skip_rows, .name_repair = "minimal") %>%
  subset(select = columns)

#rename unicode \beta to ascii "beta"
names(data)[8] <- "beta"

#extract significant and common snps:
significant_vd_snps <- data %>%
  subset(P <= p_threshold & MAF >= comman_threshold)

#write output
write.table(significant_vd_snps, output, sep = "\t", row.names = F, quote = F)
