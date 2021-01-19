#!/bin/bash -x

# Fail on errors, so-called "bash strict mode"
set -e -u -o pipefail

python3.8 -m pip install --user 'biopython == 1.77'
python3.8 -m pip install --user zstandard
python3.8 -m pip install --user tqdm

# Install tre regexp-up-to-edit-distance library
sudo apt-get update
sudo apt-get install --yes agrep libtre-dev git
if ! [ -d tre ]; then
  git clone https://github.com/ahomansikka/tre
fi
(
    cd tre/python3
    python3.8 setup.py install --user
)

python3.8 -m pip install --user sklearn
# Cython needed for compiling scikit-learn-extra
python3.8 -m pip install --user Cython
python3.8 -m pip install --user scikit-learn-extra

sudo apt-get install --yes progress
sudo apt-get install --yes htop

# the GNU `parallel` tool
sudo apt-get install --yes parallel
