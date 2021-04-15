#!/bin/bash

BASEDIR="$(dirname "$0")"

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    *)          machine="UNKNOWN"
esac

echo ${machine}

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

set -e
if [ $machine = "Mac" ]; then
	sed '/{{linux-only}}/d' "$BASEDIR"/environment.yaml  > "$BASEDIR"/environment_modified.yaml 
else
	cp "$BASEDIR"/environment.yaml "$BASEDIR"/environment_modified.yaml 
fi
#set +e
#set +eu
mamba env update -q \
	--file "$BASEDIR"/environment_modified.yaml 
#set -eu

echo CREATED the environment ${ENV}
