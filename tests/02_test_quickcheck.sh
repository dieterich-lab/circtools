#!/bin/bash

# Copyright (C) 2018 Tobias Jakobi
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

# get basic test data (for QC module)
wget https://data.dieterichlab.org/s/tzTtR6oZem5bgQd/download -O 02_quickcheck.tar.bz2
tar jxvf 02_quickcheck.tar.bz2

# change into working dir
cd 02_quickcheck/

# execute quickcheck
circtools quickcheck -d dcc_out/ -s logs/ -l AA,BB -g 1,2,1,2 -R 4,5,6,7,8,9,10,11,12,13,14,15,20,21
