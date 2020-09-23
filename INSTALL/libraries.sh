#!/bin/bash -x

# Fail on errors, so-called "bash strict mode"
set -e -u -o pipefail

python3.6 -m pip install --user 'biopython == 1.77'
python3.6 -m pip install --user zstandard
python3.6 -m pip install --user tqdm

# Install tre regexp-up-to-edit-distance library
sudo apt-get update
sudo apt-get install --yes agrep libtre-dev git
git clone https://github.com/ahomansikka/tre
(
    cd tre/python3
    python3.6 setup.py install --user
)

# the GNU `parallel` tool
sudo apt-get install --yes parallel
