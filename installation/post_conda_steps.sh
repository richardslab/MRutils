#!/bin/bash

## run post-conda steps
BASEDIR="$(dirname "$0")"

# this PAT doesn't need any permissions.
#GITHUB_PAT=$1 
echo RUNNING post-conda steps

# shellcheck source=/dev/null
set +e \
  && PS1='$$$ ' \
  && . "$(conda info --base)"/etc/profile.d/conda.sh \
  && conda activate base 
set -e

conda activate VitaminD_MR

#if [[ grep -n PAT_GITHUB ~/.Renviron ]] ; then
#	echo "GITHUB_PAT=${GITHUB_PAT}" >> ~/.Renviron
#fi

R --no-save < "$BASEDIR"/post_conda_steps.R

