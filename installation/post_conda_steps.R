#!/usr/bin/env R

Sys.setenv(TAR = 'env tar')
install.packages(c(
	"devtools", 
	"MendelianRandomization", 
	"mr.raps", 
	"rstudio"), repos = "https://cloud.r-project.org") # for some reason, now I need version 1.2.1 and that's not available via conda....

install <- function(repo, ref, dependencies = FALSE){
	devtools::install_github(repo,
		ref = ref,
		quiet = FALSE,
		dependencies = dependencies,
		upgrade = FALSE)
}

install("rondolab/MR-PRESSO", "cece763b47e59763a7916974de43c7cb93843e41")
install("WSpiller/RadialMR", "25aa8a5e64d1318c9d40646b1fd2d2847fb0a0d7")
install("gqi/MRMix", "56afdb2bc96760842405396f5d3f02e60e305039")
install("mrcieu/ieugwasr", "ff3ef11f4514ba0431bbe19635704ab1428ed129")
install("MRCIEU/MRInstruments","0.3.3")

# TODO: there are several libraries that I was unable to install via neither 
# conda nor github, so they will be installed at this point as dependencies. 
# it would be better to fix this, but I was unsuccessful

# devtools::install_version("meta",version = "4.18-0",repos = "http://cran.us.r-project.org")
# tmvnsim_1.0-2
# CompQuadForm_1.4.3
# mnormt_2.0.2
# meta_4.18-0
# psych_2.1.3
# cowplot_1.1.1
# gridExtra_2.3

install("MRCIEU/TwoSampleMR", "0.5.4", dependencies = NA)


