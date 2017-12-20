#!/usr/bin/env perl

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

# from https://www.biostars.org/p/104694/

use strict;
use warnings;

use lib "$ENV{HOME}/work/src/ensembl-api/gbioperl-1.6.1";
use lib "$ENV{HOME}/work/src/ensembl-api/gensembl/modules";
use lib "$ENV{HOME}/work/src/ensembl-api/gensembl-compara/modules";
use lib "$ENV{HOME}/work/src/ensembl-api/gensembl-variation/modules";
use lib "$ENV{HOME}/work/src/ensembl-api/gensembl-funcgen/modules";


use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

## Load the databases into the registry
$registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous'
);

## Get the gene adaptor for human
my $adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );


# Retrieve slices of every chromosome in the database
my @slices = @{ $adaptor->fetch_all() };

foreach my $gene (@slices){

    ## Get all transcripts for the gene
    my $transcripts = $gene->get_all_Transcripts;

    foreach my $transcript(@$transcripts){
#        print $transcript->display_id()."\t";
          my $introns = $transcript->get_all_Introns;
            foreach my $intron(@$introns) {

                print $intron->seq_region_name."\t";
                print "ensembl\t";
                print "intron\t";
                print $intron->start."\t";
                print $intron->end."\t";
                print ".\t";
                print $intron->strand."\t";
                print ".\t";
                print $gene->external_name."\n";


        }
    }
}

