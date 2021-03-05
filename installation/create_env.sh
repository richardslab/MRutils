#!/bin/bash

ENV=vitaminD_test
set +eu \
  && PS1=dummy \
  && . $(conda info --base)/etc/profile.d/conda.sh \
  && conda activate base \
  && conda install -y -c conda-forge mamba 

set -eu;

conda env remove ${ENV} || echo "couldn't remove environment ${ENV}"
conda create -y -n ${ENV} 
mamba env update \
	--name ${ENV} \
	--file installation/environment.yaml

conda activate ${ENV}

## run post-conda steps

R < post_conda_steps.R

