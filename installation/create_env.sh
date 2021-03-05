#!/bin/bash


ENV=vitaminD_MR
set +eu \
  && PS1=dummy \
  && conda init \
  && . $(conda info --base)/etc/profile.d/conda.sh \
  && conda activate base \
  && conda install -y -c conda-forge mamba 

set -eu;

conda env remove ${ENV} || echo "couldn't remove environment ${ENV}"
conda create -y -n ${ENV} 

set +eu
mamba env update \
	--file installation/environment.yaml
set -eu

echo CREATED the environment ${ENV}

conda activate ${ENV}

## run post-conda steps
echo RUNNING post-conda steps

R < post_conda_steps.R

