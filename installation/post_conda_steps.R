#!/usr/bin/env Rscript


#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (lock file).n", call.=FALSE)
} 
project_root <- args[1]

cat("project is", project_root,"\n")

packages_to_install <- c("dplyr","ggplot2","argparse","tidyverse","LDlinkR","devtools","openxlsx","mgsub","validate","testthat")

renv::init(project=project_root, restart=TRUE, force=TRUE)

renv::install(project=project_root, packages_to_install)

renv::install(project=project_root, "MRCIEU/TwoSampleMR")

renv::hydrate(project=project_root)

renv::snapshot(project=project_root)
