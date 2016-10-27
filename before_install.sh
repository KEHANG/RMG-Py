#!/bin/bash

# Set up anaconda
wget http://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b -p $HOME/miniconda
export PATH=$HOME/miniconda/bin:$PATH
export PYTHONPATH=$TRAVIS_BUILD_DIR/RMG-Py:$PYTHONPATH

# Update conda itself
conda update --yes conda

cd ..
git clone https://github.com/ReactionMechanismGenerator/RMG-database.git
cd RMG-Py

conda env create -f environment_linux.yml