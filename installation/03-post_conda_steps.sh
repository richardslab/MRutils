#!/bin/bash

## run post-conda steps
BASEDIR="$(dirname "$0")"

echo RUNNING post-conda steps

set +e \
  && PS1='$$$ ' \
  && . "$(conda info --base)"/etc/profile.d/conda.sh \
  && conda activate base 
#set -e

conda activate VitaminD_MR

R --no-save < "$BASEDIR"/post_conda_steps.R

