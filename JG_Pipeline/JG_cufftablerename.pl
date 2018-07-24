#!/usr/bin/perl -w

#######################################################################
#  Copyright 2015 John Garbe
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

=head1 NAME

cufftablerename.pl - Relabel sample names in a table produced by cufflinks

=head1 SYNOPSIS

cufftablerename.pl genes.fpkm_table samplenames.txt > genes.fpkm_table.renamed

=head1 DESCRIPTION

Given a *.fpkm_table file from cuffdiff2, and a file containing sample names, convert sample names in the fpkm_table header from the defulat cuffdiff q0_0 style to the names provided. The sample name file has two columns with no header. The first column is the cuff id (q1, q2, etc) and the second column is the sample name to replace the cuff id. 

=cut

use Getopt::Std;
use Pod::Usage;

our ($opt_h);

my $usage = "usage: cufftablerename.pl genes.fpkm_table samplenames.txt > genes.fpkm_table.renamed\n";
die $usage unless getopts('h');
pod2usage(q(-verbose) => 3) if ($opt_h);
die $usage unless @ARGV >= 1;

my $inputFile = $ARGV[0];
my $sampleFile = $ARGV[1];

### read in sample name file ###
open SFILE, "<", $sampleFile or die "Cannot open sample name file $sampleFile: $!\n";
while (<SFILE>) {
    chomp;
    ($group, $sample) = split;
    $replicates{$group} = 0 if (not defined($replicates{$group}));
    $name = $group . "_" . $replicates{$group};
    $label{$name} = $sample;
    $replicates{$group}++;
}

### read in fpkm_table file ###
open (INPUT, "<", $inputFile) || die("Could not open file $inputFile \n");
# header looks like this:
# tracking_id q1_0 q2_1 q3_0 ...
my $header = <INPUT>;
chomp $header;
@header = split /\t/, $header;

# read in rest of file
while (<INPUT>) {
    chomp;
    @line = split /\t/;
    $tracking_id{$line[0]} = 1;
    for $i (1..$#line) {
	$sample = $header[$i];
#	print "$i: :$sample:\n";
	$sample = $label{$sample} if %label;
	$fpkm{$sample}{$line[0]} = $line[$i];
	
#	$sample = $line[1] . "-" . $line[2];
#	$sample = $label{$sample} if %label;
#	$fpkm{$sample}{$line[0]} = $line[6];
#	$expressed{$sample}++ if ($line[6] >= $minexpression);
    }
}
close(INPUT);

@sorted_tracking_id = sort( keys(%tracking_id));

### print out data with new sample names ###
# print out header
print "$header[0]";
foreach $sample (sort( keys %fpkm)) {
    print "\t$sample";
}
print "\n";
foreach $tracking_id (@sorted_tracking_id) {
    print "$tracking_id";
    foreach $sample (sort( keys %fpkm)) {
	print "\t$fpkm{$sample}{$tracking_id}";
    }
    print "\n";
}
