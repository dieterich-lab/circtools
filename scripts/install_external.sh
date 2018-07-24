#!/bin/bash

# Copyright (C) 2017 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

function detect_os {
    unameOut="$(uname -s)"
    case "${unameOut}" in
        Linux*)     machine=Linux;;
        Darwin*)    machine=Mac;;
        CYGWIN*)    machine=Cygwin;;
        MINGW*)     machine=MinGw;;
        *)          machine="UNKNOWN:${unameOut}"
    esac
    echo ${machine}
}

function install_bedtools {
    cd /tmp/
    wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
    tar -zxvf bedtools-2.27.1.tar.gz
    cd bedtools2
    make
    mkdir -p  $HOME/.local/bin/
    cp -av bin/* $HOME/.local/bin/
    # mkdir -p  $HOME/.local/share/bedtools/
    # cp genomes -av $HOME/.local/share/bedtools/
    # rm /tmp/bedtools2 -rf
}

# install statsmodels first, does not work in setup.py due to
# https://github.com/dieterich-lab/circtools/issues/55
#pip3 install statsmodels

# install dependencies for R first
if [ "$TRAVISBUILD" ]; then
  Rscript scripts/install_R_dependencies.R
else
  Rscript scripts/install_R_dependencies.R
fi

BEDTOOLS=`which bedtools`
OS=`detect_os`

if [ "$OS" = "Mac" ]; then
  brew install libgit2
fi

if [ $BEDTOOLS ]; then

    # get current version of bedtools
    VERSION=`bedtools --version | cut -f 2 -d '.'`

    # we want to have >= 27 in order to work correctly
    if [ "$VERSION" -lt "27"  ]; then
        install_bedtools
    fi
else
     install_bedtools
fi

if [ "$OS" = "Mac" ]; then

  echo "checking for libgit2"
  if ! [[ `brew ls --versions libgit2` ]]; then
    brew install libgit2
  fi

  echo "checking for R"
  if ! [[ `brew ls --versions R` ]]; then
    brew install R
  fi

  echo "checking for python3"
  if ! [[ `brew ls --versions python@3` ]]; then
    # this is the formula for python 3.6
    # python 3.7 currently does not work with pysam
    brew install https://raw.githubusercontent.com/Homebrew/homebrew-core/f2a764ef944b1080be64bd88dca9a1d80130c558/Formula/python.rb
  fi

  echo "checking for python2"
  if ! [[ `brew ls --versions python@2` ]]; then
    brew install python@2
  fi

fi

# install DCC
cd /tmp/
git clone https://github.com/dieterich-lab/DCC.git
cd DCC
if [ "$OS" = "Mac" ]; then
  python2 setup.py install
else
  python2 setup.py install --user
fi

# install FUCHS
cd ..
git clone https://github.com/dieterich-lab/FUCHS.git
cd FUCHS

if [ "$OS" = "Mac" ]; then
  python2 setup.py install
else
  python2 setup.py install --user
fi

# remove all temporary files
#rm /tmp/FUCHS/ -rf
#rm /tmp/DCC/ -rf
