#!/bin/bash

## run post-conda steps
BASEDIR="$(dirname "$0")"


echo BASEDIR="$BASEDIR"
echo RUNNING post-conda steps

# shellcheck source=/dev/null
set +e \
  && PS1='$$$ ' \
  && . "$(conda info --base)"/etc/profile.d/conda.sh \
  && conda activate base 

set -e

pushd "${BASEDIR}/../" || exit 1
PROJECT="$(pwd)"
echo PROJECT="$PROJECT"
popd || exit 1


conda activate MRutils
Rscript "$BASEDIR"/post_conda_steps.R "$PROJECT"

# make sure that Rstudio will have the right path
echo "PATH=$PATH" > "${BASEDIR}/../.Renviron"