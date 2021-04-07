#!/bin/bash 

curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
rm Miniconda3-latest-Linux-x86_64.sh
conda init