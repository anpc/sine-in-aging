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

# pip that installs libraries specifically for python3.8
# ------------------------------------------------------
# If 3.8 is the main "system python" (in 20.04 it is):
#   + We will get it via `python3-pip` ubuntu package.
#     Has no `pip3.8` binary; OK, `python3.8 -m pip` is better anyway.
#   - CAN'T get it via `ensurepip` module
#     ("ensurepip is disabled in Debian/Ubuntu for the system python.")
# If it's not "system python":
#   ? `python3-pip` package MAY OR MAY NOT install it for 3.8.
#   ? `python3.8-venv` package MAY OR MAY NOT install pip;
#   + BUT `python3.8-venv` contains `ensurepip` module that can install pip.
# Both ways may give an old pip, but whatever, it works.

sudo apt-get install --yes python3-pip python3.8-venv
if ! python3.8 -m pip --version
then
    sudo python3.8 -m ensurepip
fi

pip3.8 --version || echo "That's OK, calling 'python3.8 -m pip' is more reliable anyway."
python3.8 -m pip --version