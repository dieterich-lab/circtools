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

# get basic test data (humane genome)
wget https://data.dieterichlab.org/s/eikQFHKFstSgbrp/download -O 00_base.tar.bz2
tar jxvf 00_base.tar.bz2

# get test data for DCC
wget https://data.dieterichlab.org/s/pn7QHoQJmtD44Fo/download -O 01_dcc.tar.bz2
tar jxvf 01_dcc.tar.bz2

# get basic test data (humane genome)
wget https://data.dieterichlab.org/s/emNDzztToQoyerz/download -O chr1.gtf.bz2
bunzip2 chr1.gtf.bz2

# change into working dir
cd 01_dcc/

# execute DCC
circtools detect @samplesheet -ss -T 2 -D -an ../chr1.gtf -A ../00_base/GRCh38_85.fa -R ../00_base/GRCh38_85_repeatmasker.gtf -B @bam_files.txt -M -Nr 2 2 -fg -G -t /tmp/ -F -L 20 -k -O output

if [ "$OS" = "Mac" ]; then
  md5 output/*
else
  md5sum output/*
fi
