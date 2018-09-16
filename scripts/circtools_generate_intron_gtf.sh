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
if [ ! $# == 2 ]; then
  echo "Usage: $0 [Input GTF] [target dir e.g. /tmp/]"
  exit
fi

BEDTOOLS=`which bedtools`

if [ $BEDTOOLS ]; then
	awk '{if ($3=="gene"){print $0}}' $1 > $2/genes.only.gtf
	awk '{if ($3=="exon"){print $0}}' $1 > $2/exons.only.gtf
	bedtools subtract -a $2/genes.only.gtf -b $2/exons.only.gtf > $2/introns.gtf
	sed -i 's/\Wgene_id.*ene_name \"/\t/g' $2/introns.gtf
	sed -i 's/\"; gene_source.*//g' $2/introns.gtf
	awk 'BEGIN{OFS="\t"} {$9="gene_name \""$9"\";"; print $0}' $2/introns.gtf | uniq > $2/introns.names_fixed.gtf
	sed -i 's/\Wgene\W/\tintron\t/g' $2/introns.names_fixed.gtf
	cat $2/genes.only.gtf > $2/genes_and_introns.gtf
	cat $2/introns.names_fixed.gtf >> $2/genes_and_introns.gtf
	rm $2/introns.gtf
	rm $2/genes.only.gtf
	rm $2/exons.only.gtf
	rm $2/introns.names_fixed.gtf		
else
	echo "This script requires bedtools to be installed."
fi

