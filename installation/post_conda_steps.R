#!/usr/bin/env R

Sys.setenv(TAR = '/bin/tar')
install.packages("devtools", repos = "https://cloud.r-project.org") # for some reason, now I need version 1.2.1 and that's not available via conda....

install <- function(repo, ref){
	devtools::install_github(repo,
		ref = ref,
		quiet = FALSE,
		dependencies = FALSE,
		upgrade = "never")
}

install("rondolab/MR-PRESSO", "cece763b47e59763a7916974de43c7cb93843e41")
install("WSpiller/RadialMR", "25aa8a5e64d1318c9d40646b1fd2d2847fb0a0d7")
install("gqi/MRMix", "56afdb2bc96760842405396f5d3f02e60e305039")
install("mrcieu/ieugwasr", "ff3ef11f4514ba0431bbe19635704ab1428ed129")
install("MRCIEU/TwoSampleMR", "0.5.4")


