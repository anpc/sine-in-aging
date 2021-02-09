#!/bin/bash -x

# Fail on errors, so-called "bash strict mode"
set -e -u -o pipefail

cat /etc/lsb-release  # Just to see which Ubuntu version we're on

sudo apt-get update
sudo apt-get install --yes curl
sudo apt-get install --yes software-properties-common # for add-apt-repository

# For python 3.8 on ubuntu 16.04
# https://askubuntu.com/questions/865554/how-do-i-install-python-3-6-using-apt-get
if ! which python3.8
then
    sudo add-apt-repository --yes ppa:deadsnakes/ppa
    sudo apt-get update
fi

sudo apt-get install --yes python3.8 python3.8-dev