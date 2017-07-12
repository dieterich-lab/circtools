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

# install primer-design R part
Rscript R/install.R

# install CircTest
Rscript install_circtest.R

# install dependencies for alternative exon usage
Rscript install_depenencies_exon_usage.R

# install DCC
cd /tmp/
git clone https://github.com/dieterich-lab/DCC.git
cd DCC
python setup.py install --user

# install FUCHS
cd ..
git clone https://github.com/dieterich-lab/FUCHS.git
cd FUCHS
python setup.py install --user

# remove all temporary files
rm /tmp/FUCHS/ -rf
rm /tmp/DCC/ -rf

