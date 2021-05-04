#!/bin/bash

BASEDIR="$(dirname "$0")"

ENV=MRutils
# shellcheck source=/dev/null
set +eu \
  && PS1='$$$ ' \
  && . "$(conda info --base)"/etc/profile.d/conda.sh \
  && conda activate base \
  && conda install -y -c conda-forge mamba 

conda deactivate || echo "No active environment"
conda env remove -n ${ENV} || echo "Couldn't remove environment ${ENV}"
conda create -y -n ${ENV}  || echo "It seem that environment ${ENV} is already present"

set -ex

mamba env update -n ${ENV} -q \
	--file "$BASEDIR"/environment_modified.yaml 

echo CREATED the environment ${ENV}
