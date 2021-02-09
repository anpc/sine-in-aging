#!/bin/bash -x

# Fail on errors, so-called "bash strict mode"
set -e -u -o pipefail

cat /etc/lsb-release  # Just to see which Ubuntu version we're on

cd "$(dirname "$0")"  # the directory of this script
./python.sh
./pip.sh
./libraries.sh
