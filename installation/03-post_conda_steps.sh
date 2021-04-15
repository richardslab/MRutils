#!/bin/bash

## run post-conda steps
BASEDIR="$(dirname "$0")"

echo RUNNING post-conda steps

# shellcheck source=/dev/null
set +e \
  && PS1='$$$ ' \
  && . "$(conda info --base)"/etc/profile.d/conda.sh \
  && conda activate base 

conda activate VitaminD_MR
echo GH=$GITHUB_PAT
R --no-save < "$BASEDIR"/post_conda_steps.R

