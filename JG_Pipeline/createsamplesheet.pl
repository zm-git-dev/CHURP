#!/usr/bin/perl -w

#######################################################################
#  Copyright 2014 John Garbe
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

createsamplesheet.pl - Automatically create a samplesheet for a folder of fastq files

=head1 SYNOPSIS

createsamplesheet.pl -f folder

=head1 DESCRIPTION

Automatically generate a samplesheet suitable for use with various riss_util scripts and pipelines. The samplesheet is a Qiime mapping file if Qiime mode is enabled with the "-q" option.

 -f folder : A folder containing fastq files to process
 -o file : Name of the output samplesheet file
 -z : fastq files are gzip compressed (filenames end with fastq.gz)
 -h : Print usage instructions and exit
 -v : Print more information while running (verbose)
 -i : Assume Illumina style fastq file names
 -q : Enable Qiime mode: generate a Qiime-compatible mapping file
 -p file : A fasta file containing possible primers to look for (Qiime-mode)
 -e : Do not attempt to identify primers in the reads (Qiime-mode: for use with EMP protocol datasets, or any other 16s protocol where primers are not present in the reads)
 -s integer : Only use this many sequences from the beginning of each fastq file for primer detection (Qiime-mode, default: 4000)

=cut

############################# Main  ###############################

my $starttime = time();
my $checkpointtime = $starttime;

use Getopt::Std;
use Cwd 'abs_path';
use File::Basename;
use Pod::Usage;
#use feature "switch";
use Scalar::Util qw(looks_like_number);
#use FindBin qw($RealBin);

my %samples;
my %args;
my %stats;

our ($opt_f, $opt_o, $opt_h, $opt_p, $opt_s, $opt_v, $opt_e, $opt_q, $opt_z, $opt_i);

$usage = "USAGE: createsamplesheet.pl [-f /fastq/folder]\n";
die $usage unless getopts('f:o:hp:s:veqzi');
die $usage if ($#ARGV >= 0);
pod2usage(q(-verbose) => 3) if ($opt_h);

### parameters
$fastqfolder = abs_path($opt_f // $ENV{PWD});
die "Fastq folder $fastqfolder not found\n" if (!-e $fastqfolder);
$opt_o = abs_path($opt_o) if ($opt_o);
$outputfile = $opt_o // "samplesheet.txt";
print "Creating samplesheet file: $outputfile\n";
$verbose = $opt_v // 0; # print out more stuff
$gz = $opt_z // 0;
$qiime = $opt_q // 0;
$primerfile = $opt_p // "/home/umii/public/gopher-pipelines/1.5/resources/16s-primers.fa";
$subsample = $opt_s // 4000;
$noprimers = $opt_e // 0; # don't detect primers (emp protocol)

### process fastq directory
# do a lot of work to get short sample names from the fastq filenames
if ($gz) {
    @files = <$fastqfolder/*.fastq.gz>;
} else {
    @files = <$fastqfolder/*.fastq>;
    if ($#files < 0) {
	print "No fastq files found, looking for compressed files...\n";
	@files = <$fastqfolder/*.fastq.gz>;		
    }
    $gz = 1;
}
die "No fastq files found in $fastqfolder\n" if ($#files < 0);

foreach $file (@files) {
    my ($fname, $path) = fileparse($file);
    if ($fname =~ /_(R[12])[._]/) {
	push @fastqfiles, $file;
    } else {
	print "Ignoring file $file, unknown filename format\n";
    }
}

# Split the filenames based on illumina format
if ($opt_i) {
    %samples = ();
    foreach $file (@fastqfiles) {
	my ($fname, $path, $mysuffix) = fileparse($file, (".fastq", ".fastq.gz"));

	if ($fname =~ /(.*)_([AGCT]{3,8})(-[AGCT]{3,8})?_L00([1-8])_(R[12])_[0-9]{3}$/) {
	    # standard old-style filenames
	    print STDERR "Filename format: old-style\n" if ($verbose);
	    $sample = $1;
	    $read = $5;
	    $lane = $4;
	    $sample = "$sample.L00$lane";
	} elsif ($fname =~ /(.*)_S[\d]{1,4}_L00([1-8])_(R[12])_[0-9]{3}$/) {
	    # standard new-style filenames
	    print STDERR "Filename format: new-style\n" if ($verbose);
	    $sample = $1;
	    $read = $3;
	} elsif ($fname =~ /(.*)_S[\d]{1,4}_(R[12])_[0-9]{3}$/) {
	    # standard new-style concatenated (no lane number) filenames
	    print STDERR "Filename format: new-style concatenated\n" if ($verbose);
	    $sample = $1;
	    $read = $2;
	} elsif ($fname =~ /lane[1-8]_Undetermined_L00([1-8)_(R[12])_[0-9]{3}$/) {
	    # old-style undetermined filenames
	    print STDERR "Filename format: old-style undetermined\n" if ($verbose);
	    $sample = "lane${1}-Undetermined";
	    $read = $2;
	} elsif ($fname =~ /Undetermined_S0_(R[12])_[0-9]{3}$/) {
	    # new-style undetermined concatenated filenames
	    print STDERR "Filename format: new-style undetermined concatenated\n" if ($verbose);
	    $sample = "Undetermined";
	    $read = $1;
	} else {
	    print "Unable to parse $fname, skipping\n";
	    next;
	}

	$samples{$sample}{$read}{fastq} = $file;
	$samples{$sample}{$read}{fastqname} = $fname . $mysuffix;
    }

} else {
    # Split the filenames at the first, second, third, etc underscore character, until the split filenames are unique for each sample
    for $i (0..7) {
	%samples = ();
	$pe = 0;
	foreach $file (@fastqfiles) {
	    my ($fname, $path) = fileparse($file);
	    @parts = split /_/, $fname;
	    if ($#parts <= $i) {
		$sample = $fname;
	    } else {
		$sample = join "_", @parts[0..$i];
	    }
#	    $sample =~ s/_/-/g; # convert _ to -
	    if ($fname =~ /_(R[12])[._]/) {
		$read = $1;
		$pe = 1 if ($read eq "R2");
		$samples{$sample}{$read}{fastq} = $file;
		$samples{$sample}{$read}{fastqname} = $fname;
		
	    }
	}
	$samplenumber = keys %samples;
#       print "i: $i, $samplenumber samples; $#fastqfiles fastqfiles\n";
	if ($pe) {
	    last if (keys %samples eq (($#fastqfiles + 1) / 2));
	} else {
	    last if (keys %samples eq ($#fastqfiles + 1));
	}
    }
}

# determine number of R1 and R2 files
$r1count = 0;
$r2count = 0;
foreach $sample (keys %samples) {
    $r1count++ if ($samples{$sample}{R1});
    $r2count++ if ($samples{$sample}{R2});
}
    
if ($r2count == 0) {
    print "Single-end read dataset\n";
    $pe = 0;
} elsif ($r2count == $r1count) {
    print "Paired-end read dataset\n";
    $pe = 1;
} else {
    print "Unequal number of R1 and R2 fastq files found: $#r1count R1 fastq files, $#r2count R2 files. Mixing paired-end and single-end datasets is not supported\n";
    exit;
}

$samplenumber = keys %samples;
print "$samplenumber samples found in folder $fastqfolder\n";
die "No fastq files with BMGC format filenames could be found\n" if ($samplenumber == 0);

### Qiime only sections ###
if ($qiime) {
### Read in primer file ###
    unless ($noprimers) {
	open PFILE, "$primerfile" or die "Cannot open primer file $primerfile\n";
	$first = 1;
	$primers{""} = "undef";
	while ($line = <PFILE>) {
	    $line =~ s/\r\n?/\n/g; # get rid of stupid windows newlines
	    chomp $line;
	    
	    if ($line =~ ">") {
		$line =~ s/^>//;
		$primers{$seq} = $id unless ($first);
		push @primers, $seq unless ($first);
		$first = 0;
		$id = $line;
		$seq = "";
	    } else {
		# translate IUPAC codes to Ns - cutadapt 1.7 now supports IUPAC, yay
#	$line =~ tr/acgtACGT/N/c; 
		$seq = "$seq$line";
	    }
	}
	$primers{$seq} = $id unless ($first);
	push @primers, $seq unless ($first);
	$primernumber = (keys %primers) - 1; # Don't count the "undefined" primer
	die "No valid primers found in file $primerfile\n" if ($primernumber < 1);
	print "$primernumber primers found in file $primerfile\n";
    }

    ### Identify Primers ###
    unless ($noprimers) {
	print "\nSample\tPrimer\tFraction";
	print "\tPrimer2\tFraction" if ($pe);
	print "\n";
    }
    foreach $sample (keys %samples) {
#	print "$sample, $samples{$sample}{R1}, $samples{$sample}{R2}\n";
	if ($noprimers) {
	    $stats{$sample}{LinkerPrimerSequence} = "";
	    $stats{$sample}{primernamer1} = "noprimer";
	    if ($pe) {
		$stats{$sample}{ReversePrimer} = "";
		$stats{$sample}{primernamer2} = "noprimer";
	    }
	} else {
	    
	    my $forwardprimer = "";
	    my $reverseprimer = "";
	    
	    print "$sample";
	    ($forwardprimer, $pct1) = &primerdetect($samples{$sample}{R1}{fastq});
	    print "\t$primers{$forwardprimer}\t$pct1";
	    $stats{$sample}{LinkerPrimerSequence} = $forwardprimer;
	    $stats{$sample}{primernamer1} = $primers{$forwardprimer};
	    if ($pe) {
		($reverseprimer, $pct2) = &primerdetect($samples{$sample}{R2}{fastq});
		print "\t$primers{$reverseprimer}\t$pct2";
		$stats{$sample}{ReversePrimer} = $reverseprimer;
		$stats{$sample}{primernamer2} = $primers{$reverseprimer};
	    }
	    print "\n";
	}
    }
}

# Add additional attributes
foreach $sample (sort keys %samples) {
    $stats{$sample}{Group} = substr $sample, 0,2;
    $groups{$stats{$sample}{Group}} = 1; # keep track of groups
    $stats{$sample}{fastqR1} = $samples{$sample}{R1}{fastqname};
    $stats{$sample}{fastqR2} = $samples{$sample}{R2}{fastqname} if ($pe);
    $stats{$sample}{Description} = "Sample-$sample";
    $sanitized = $sample;
    $sanitized =~ s/[^a-z0-9]/\./ig;
    $stats{$sample}{"#SampleID"} = $sanitized;
    $stats{$sample}{BarcodeSequence} = "";
}

# If there are not at least two groups, regenerate the group attribute so that there are two groups
if (keys %groups < 2) {
    $count = 0;
    foreach $sample (sort keys %samples) {
	$count++;
	if ($count % 2 == 0) {
	    $stats{$sample}{Group} = "A";
	} else {
	    $stats{$sample}{Group} = "B";
	}
    }
}
    

### Print samplesheet ###
# Determine columns to print
if ($qiime) {
    if ($pe) {
	@columns = ("#SampleID", "BarcodeSequence", "LinkerPrimerSequence", "ReversePrimer", "fastqR1", "fastqR2", "Group", "Description");
    } else {
	@columns = ("#SampleID", "BarcodeSequence", "LinkerPrimerSequence", "fastqR1", "Group", "Description");
    }
} else {
    if ($pe) {
	@columns = ("#SampleID", "fastqR1", "fastqR2", "Group", "Description");
    } else {
	@columns = ("#SampleID", "fastqR1", "Group", "Description");
    }
}

open OFILE, ">$outputfile" or die "Cannot open file $outputfile for writing: $!\n";

foreach $column (@columns) {
    print OFILE "$column\t";
}
print OFILE "\n";

foreach $sample (sort keys %stats) {
    foreach $column (@columns) {
	if (defined ($stats{$sample}{$column})) {
	    print OFILE $stats{$sample}{$column} . "\t";
	} else {
	    print OFILE "undef\t";
	}
    }
    print OFILE "\n";
}
close OFILE;

### Validate samplesheet ###
if (0) { # turned this off so you can use it for riss_util scripts without loading qiime
    if ($noprimers) {
	$result = `validate_mapping_file.py -m $outputfile --not_barcoded --disable_primer_check`;
    } else {
	$result = `validate_mapping_file.py -m $outputfile --not_barcoded`;
    }
    print "$result\n";
}

#&runtime;

exit;

################################ Helper subs ###############################

# determine the primer sequence present at the beginning of the reads in 
# in a fastq file
sub primerdetect {
    my ($file) = @_;

    # build list of primers
    my $primers = "";
    foreach $primer (@primers) {
	next if ($primer eq '');
	$primers .= "-g $primer ";
    }

#    print "module load python; module load cutadapt; head -n $subsample $file | cutadapt - $primers -O 10 -e 0.0 --match-read-wildcards 2>&1 1>/dev/null\n";

    if ($file =~ /\.gz$/) {
	@results =  `module load python; module load cutadapt/1.7.1; gunzip -c $file | head -n $subsample | cutadapt - $primers -O 10 -e 0.0 --match-read-wildcards 2>&1 1>/dev/null\n`;
    } else {
	@results =  `module load python; module load cutadapt/1.7.1; head -n $subsample $file | cutadapt - $primers -O 10 -e 0.0 --match-read-wildcards 2>&1 1>/dev/null\n`;
    }
#    @results =  `module load python; module load cutadapt; head -n $subsample $file | cutadapt - $primers -O 10 -e 0.0 --match-read-wildcards 2>&1 1>/dev/null\n`;

    print "Scanning for primers in file $file\n" if ($verbose);
    $topseq = "";
    $toppercent = 0;
    foreach $line (@results) {
	if ($line =~ /Processed reads:\s+(\d+)/) {
	    $numreads = $1;
	    return ('', 0) if ($numreads == 0); # avoid divide by zero from empty files
	}
	if ($line =~ /^Sequence: (\w+);/) {
	    $seq = $1;
	    $line =~ /; Trimmed: (\d+) times./;
	    $count = $1;
	    $percent = $count / $numreads;
	    print "Primer: $primers{$seq}, Sequence: $seq, Count: $count, Fraction: $percent\n" if ($verbose);
	    if ($percent > $toppercent) {
		$toppercent = $percent;
		$topseq = $seq;
	    }
	}
    }
    
#    print "$topseq\t$toppercent\n";
    if ($toppercent < 0.001) {
	print "No primers identified in file $file\n";
	print "Best primer: $primers{$topseq}, $topseq, fraction: $toppercent\n";
	return ("", 0);
    }
    return ($topseq, $toppercent);
}

# this truncates instead of rounds, fix it
sub round10 {
    return int($_[0] * 10) / 10;
}

# this truncates instead of rounds, fix it
sub round100 {
    return int($_[0] * 100) / 100;
}

# print how long we've been running
sub runtime {
    $now = time();
    my $runtime = $now - $starttime;
    my $checktime = $now - $checkpointtime;
    print "Checkpoint time: " . &formattime($checktime) . "\tTotal time: " . &formattime($runtime) . "\n";
    $checkpointtime = $now;
}

# print out checkpoint time info
sub formattime {
    my ($time) = @_;
    if ($time < 60) {
	return "$time seconds";
    } elsif ($time < 7200) { # 2 hours
	my $minutes = &round10($time / 60);
	return "$minutes minutes";
    } else {
	$hours = &round10($time / 3600);
	return "$hours hours";
    }
}
