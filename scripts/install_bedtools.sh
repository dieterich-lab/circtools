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

# this is a workaround in order to be able to install R packages as regular user

# This script downloads the installation zip of OligoArrayAux and installs it with --prefix=$HOME/.local/
# This is consistent with the circtools installation

cd /tmp/
wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
tar -zxvf bedtools-2.26.0.tar.gz
cd bedtools2
make
mkdir -p  $HOME/.local/bin/
cp bin/* $HOME/.local/bin/
mkdir -p  $HOME/.local/share/bedtools/
cp genomes -av $HOME/.local/share/bedtools/