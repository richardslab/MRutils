# MR utils 
Yossi's attempt to encapsulate the environment and functions needed to prepare data for, and use the R library TwoSampleMR

Can be used as a starting point for a Mendelian Randomization analysis.
----
# What is provided here?
This repository contains three main components:

1. Scripts for creating a reproducible computing environment (an R installation, the TwoSampleMR library, etc.) in which one can use R to perform Mendelian Randomization analyses (e.g. on your machine, but also provided in a [docker image](https://hub.docker.com/repository/registry-1.docker.io/richardslab/mr_util)
2. R functions for performing MR-related activities (getting proxies, resolving rsids, etc) currently contained in the `utils.R` file 

# How can I use this?
There are three options:
1. Install the R library TwoSampleMR and some others on your own, and simply use the methods in utils.R 
2. Use the scripts in the [installation] directory (they are numbered, run them in order)
3. Run inside the  docker container `richardslab/mr_utils`
```bash
docker run --rm -it -v $(pwd):/work richardslab/mr_utils:latest Rscript /work/relative/path/to/r/file -arg1 value1 -arg2 value2
```

This is work in progress. please visit often.
