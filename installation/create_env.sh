#!/bin/bash

BASEDIR=$(dirname $0)

ENV=vitaminD_MR

set +eu \
  && PS1=dummy \
  && . $(conda info --base)/etc/profile.d/conda.sh \
  && conda activate base \
  && conda install -y -c conda-forge mamba 

set -eu;

conda deactivate || echo "no active environment"
conda remove -n ${ENV} || echo "couldn't remove environment ${ENV}"
conda create -y -n ${ENV}  || echo "it seem that environment ${ENV} is already present"


mamba env update \
	--file $BASEDIR/environment.yaml


echo CREATED the environment ${ENV}

