#!/bin/bash --login
set -eu

# activate conda environment and let the following process take over
conda activate VitaminD_MR
exec "$@"