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

function install_bedtools {
    cd /tmp/

    if [[ "$OSTYPE" == "linux-gnu" ]]; then
        wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        curl -O https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz
    else
        echo "Sorry, this OS type not supported. Please contact circtools@dieterichlab.org for help."
    fi

    tar -zxvf bedtools-2.27.1.tar.gz
    cd bedtools2
    make
    mkdir -p  $HOME/.local/bin/
    cp -av bin/* $HOME/.local/bin/
    # mkdir -p  $HOME/.local/share/bedtools/
    # cp genomes -av $HOME/.local/share/bedtools/
    rm /tmp/bedtools2 -rf
}


# install dependencies for R first
Rscript scripts/install_R_dependencies.R

BEDTOOLS=`which bedtools`

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

# install DCC
cd /tmp/
git clone https://github.com/dieterich-lab/DCC.git -b devel
cd DCC
python3 setup.py install

# install FUCHS
cd ..
git clone https://github.com/dieterich-lab/FUCHS.git -b FUCHS-0.3
cd FUCHS
python3 setup.py install


# remove all temporary files
rm /tmp/FUCHS/ -rf
rm /tmp/DCC/ -rf

