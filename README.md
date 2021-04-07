# VitaminD_COVID_MR_reproduce
Yossi's attempt to repoduce Guilloumme's vitamin_d_to_covid_outcome MR paper

Can be used as a starting point for a Mendelian Randomization analysis.
----
# What is provided here?
This repository contains three main components:

1. Scripts for creaeing a reproducible computing environment (an R installation, the TwoSampleMR library, etc.) in which one can use R to perform Mendelian Randomization analyses (either on your machine, but also provided in a [docker image](https://hub.docker.com/repository/registry-1.docker.io/richardslab/vitamin_d_mr_test)
2. R functions for performing MR-related activities (getting proxies, resolving rsids, etc) currently contained in the `utils.R` file 
3. An example MR analysis that re-does the MR part of Butler-Laporte's [Vitamin-D ==> COVID-19 outcome MR paper](https://www.medrxiv.org/content/10.1101/2020.09.08.20190975v5)

# How can I use this?
In order to use this package you can pull the docker image and run your analysis against that image. For example, to run R interactive environment from within the environment provided by the docker image you would to something like this:

```bash
docker run --rm -it -v $(pwd):/work richardslab/vitamin_d_mr_test:latest R
```
...
