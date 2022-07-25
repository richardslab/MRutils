#!/usr/bin/env Rscript



args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (lock file).n", call.=FALSE)
} 
project_root <- args[1]

cat("project is", project_root,"\n")

packages_to_install <- c("dplyr","ggplot2","argparse","tidyverse","LDlinkR","openxlsx","mgsub","validate", "testthat")
packages_needed_for_dev <- c("devtools")

cat("R-ENV init\n")
  renv::init(project=project_root, restart=TRUE, force=TRUE)

cat("R-ENV install packages\n")
for (package in packages_to_install){
  cat("R-ENV installing: ", package, "\n")
  renv::install(project=project_root, package)
}

cat("R-ENV install github TwoSampleMr\n")
renv::install(project=project_root, "MRCIEU/TwoSampleMR")

# cat("RNV hydrate")
# renv::hydrate(project=project_root)

# cat("R-ENV snapshot\n")
# renv::snapshot(project=project_root, directory="R")
