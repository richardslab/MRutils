#!/bin/bash

ENV=vitaminD_covid_MR2
set +eu \
  && PS1=dummy \
  && . $(conda info --base)/etc/profile.d/conda.sh \
  && conda activate base;
set -eu;


mamba create --yes \
	--name ${ENV} \
	--channel bioconda \
	--channel conda-forge \
	--channel r \
	--file requirements.txt

conda activate ${ENV}

R < post_conda_steps.R

