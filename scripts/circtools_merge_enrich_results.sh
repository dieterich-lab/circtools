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


# check if we have 2 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [Input path to enrich results] [Number of iterations (for regex)] [target dir e.g. /tmp/] [file name]"
  exit
fi

awk -v OFS='\t' -F '\t' '{print FILENAME,$0}' $1/*.csv |  sed "s/_${2}_.*.csv//g" | sed "s~$1~~" | sed "s~/~~" | grep -v circRNA_host_gene > $3/$4
