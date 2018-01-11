#!/usr/bin/env perl

###
#
# Script takes CircSkipJunctions from DCC folder and
# write a BED file of linear splice junctions that intersect with DCC predicted circles
#
###

use strict;
use warnings;
my $CircSkipJunctions = shift @ARGV;

open(IN, $CircSkipJunctions) || die "Could not open CircSkipJunctions";
my $line = <IN>;
chomp $line;
print 'track name="CircSkip_SJ" description="New Linear Splice Junction after BS event" color="red"\n';

my %hash;
my %hashName;

# aggregate by SJ coordinates and not by Circle !

while ($line = <IN>) {
    my @data = split(/\t+/, $line);

    my $CircID = join("_", @data[0 .. 2]);

    my @skips = grep {/\:/} @data[3 .. $#data];
    my $sampleID = 0;

    while (my $record = shift @skips) {
        #split if multiple entries by ";"

        map {/(\S+)\:(\d+)/;
            $hash{$1}{$sampleID} = $2;
            $hashName{$1}{$CircID} = 1;} split(/\;/, $record);
        $sampleID++;
    }
}
close(IN);

foreach my $key (keys %hash) {
    if ($key) {

        $key =~ /(\S+)\:(\d+)\-(\d+)/;
        my $SumOfSJ = 0;

        map {$SumOfSJ += $_} values %{$hash{$key}};
        print $1, "\t", $2, "\t", $3, "\tCirc_", join("|", keys %{$hashName{$key}}), "\t", $SumOfSJ, "\t.\n";

    }
}
