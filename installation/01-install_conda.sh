#!/bin/bash 

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    *)          machine="UNKNOWN"
esac

if [ "$machine" = "Mac" ]; then
	curl -LO https://repo.anaconda.com/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
	bash  ./Miniconda2-latest-MacOSX-x86_64.sh -b -f -p "${HOME}/miniconda"
	rm Miniconda2-latest-MacOSX-x86_64.sh
	PATH="${HOME}/miniconda/bin:${PATH}"
	conda init
fi



if [ "$machine" = "Linux" ]; then
	set -xe
	curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash ./Miniconda3-latest-Linux-x86_64.sh -b -f -p ${HOME}/miniconda 
	rm Miniconda3-latest-Linux-x86_64.sh
	PATH="/miniconda/bin:${PATH}"
	conda init
fi