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
wget "http://unafold.rna.albany.edu/cgi-bin/OligoArrayAux-download.cgi?oligoarrayaux-3.8.tar.bz2" -O oaa.tar.bz2
tar xvf oaa.tar.bz2
cd oligoarrayaux-3.8/
./configure --prefix=$HOME/.local/
make
make install
rm /tmp/oligoarrayaux-3.8/ -rf