#!/bin/bash -x

# Fail on errors, so-called "bash strict mode"
set -e -u -o pipefail

# pip that installs libraries specifically for python3.8
# ------------------------------------------------------
#
# Previous theory for getting pip from apt packages:
#
# If 3.8 is the main "system python" (in 20.04 it is):
#   + We will get it via `python3-pip` ubuntu package.
#     Has no `pip3.8` binary; OK, `python3.8 -m pip` is better anyway.
#   - CAN'T get it via `ensurepip` module
#     ("ensurepip is disabled in Debian/Ubuntu for the system python.")
# If it's not "system python":
#   ? `python3-pip` package MAY OR MAY NOT install it for 3.8.
#   ? `python3.8-venv` package MAY OR MAY NOT install pip;
#   + BUT `python3.8-venv` contains `ensurepip` module that can install pip.
#
# Both ways may give an old pip, and unfortunately that's become a problem -
# at least 16.04 gives a pip so old that it no longer works in 3.8.
# 
# => So using https://stackoverflow.com/a/60046826 advice for installing
#    a per-python pip.

sudo apt-get install --yes python3.8-venv python3-setuptools
# `pip list` is a stronger check that it's actually working.
if ! (python3.8 -m pip --version && python3.8 -m pip list)
then
    #sudo python3.8 -m ensurepip
    sudo python3.8 -m easy_install pip
fi

python3.8 -m pip --version
python3.8 -m pip list
