#!/usr/bin/env R

Sys.setenv(TAR = '/bin/tar')
install.packages("devtools",repos="https://cloud.r-project.org") # for some reason, now I need version 1.2.1 and that's not available via conda....
devtools::install_github("MRCIEU/TwoSampleMR", 
                         ref = "0.5.4", 
                         quiet = TRUE, 
                         upgrade = "always")

