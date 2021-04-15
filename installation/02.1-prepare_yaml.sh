#!/bin/bash

BASEDIR="$(dirname "$0")"

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    *)          machine="UNKNOWN"
esac

echo ${machine}

if [ $machine = "Mac" ]; then
  sed '/{{linux-only}}/d' "$BASEDIR"/environment.yaml  > "$BASEDIR"/environment_modified.yaml 
else
  cp "$BASEDIR"/environment.yaml "$BASEDIR"/environment_modified.yaml 
fi
