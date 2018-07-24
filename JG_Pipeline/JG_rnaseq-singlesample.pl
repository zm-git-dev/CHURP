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

rnaseq-singlesample.pl - Analyze one RNA_Seq sample, from raw fastq file(s) to abundance estimation, using Hisat2 and Cufflinks

=head1 SYNOPSIS

USAGE: rnaseq-singlesample.pl --hisat2index index --gtffile file [--outputfolder folder] --R1 file [--R2 file]

=head1 DESCRIPTION

Analyze one RNA-Seq sample:

 --hisat2index index : A hisat2 index basename (for running hisat2)
 --bowtie2index index : A bowtie2 index basename (for running tophat2)
 --gtffile string : GTF annotation file
 --stranded : TruSeq Stranded dataset
 --maskfile : GTF file for Cuffquant's --mask-file option
 --outputfolder folder : The output folder
 --threads integer : Number of processors to use (number of threads to run)
 --verbose : Print more information while running
 --subsample integer : Subsample the specified number of reads from each sample. 0 = no subsampling (default = 0)
 --qualitycontrol : Enable quality-control: use trimmomatic to trim adapter sequences, trim low quality bases from the ends of reads, and remove short sequences
 --headcrop integer : crop the specified number of bases off the front of each read, useful for picogreen libraries
 --extraoptionsfile file : File with extra options for trimmomatic, tophat, cuffquant, or featurecounts
 --help : Print usage instructions and exit
=cut

##################### Main ###############################

use FindBin;                  # locate this script directory
use lib "$FindBin::RealBin";  # look for modules in the script directory
use Pipelinev4;               # load the pipeline module

use Getopt::Long;
use Pod::Usage;
use Cwd 'abs_path';

my %args;
my %stats;

### Initialize ###
&init(\%args, \%stats);

### Subsample ###
runtime;
subsample(\%args, \%stats);

### QC ###
#runtime;
#qc(\%args, \%stats);

### FastQC ###
runtime;
fastqc(\%args, \%stats);

if ($args{qualitycontrol}) {
    ### Trimmomatic ###
    runtime;
    trimmomatic(\%args, \%stats);

    ### FastQC ###
    runtime;
    fastqc(\%args, \%stats, "fastqc-trim");
}

if ($args{hisat2index}) {
    ### Hisat2 ###
    runtime;
    hisat2(\%args, \%stats);
} else {
    ### Tophat2 ###
    runtime;
    tophat2(\%args, \%stats);
}

### Picard cleansam ###
#runtime;
#&cleansam(\%args, \%stats);

### Samtools sort ###
runtime;
samtoolssortandindexbam(\%args, \%stats);

### Picard insert metrics ###
runtime;
insertsize(\%args, \%stats);

### Subread featureCounts ###
runtime;
featurecounts(\%args, \%stats);

### Cuffquant ###
runtime;
cuffquant(\%args, \%stats);

### Write stats to file ###
runtime;
writestats(\%args, \%stats);

`touch $args{scratchfolder}/Complete`;

print "Finished\n";
runtime;

exit;

##################### Initialize ###############################
sub init {
    print "Running Sample Initialize\n";

    my ($args, $stats) = @_;

    ### Usage ###
    $usage = "USAGE: rnaseq-singlesample.pl --hisat2index index --gtffile file [--outputfolder folder] --R1 file [--R2 file]\n";

    # set defaults - add additional parameter defaults here
    $args->{threads} = $ENV{"PBS_NUM_PPN"} // 1;
    $args->{scratchfolder} = $ENV{PWD} . "/rnaseq2-singlesample";
    GetOptions("help" => \$args->{help},
               "verbose" => \$args->{verbose},
	       "R1=s" => \$args->{R1},
	       "R2=s" => \$args->{R2},
               "threads=i" => \$args->{threads},
               "outputfolder=s" => \$args->{scratchfolder},
               "extraoptionsfile=s" => \$args->{extraoptionsfile},
               "subsample=i" => \$args->{subsample},
               "qualitycontrol" => \$args->{qualitycontrol},
	       "hisat2index=s" => \$args->{hisat2index},
               "bowtie2index=s" => \$args->{bowtie2index},
	       "gtffile=s" => \$args->{gtffile},
	       "maskfile=s" => \$args->{maskfile},
               "stranded" => \$args->{stranded},	       
	       "headcrop=i" => \$args->{headcrop},
        ) or die $usage;
    pod2usage(q(-verbose) => 3) if ($args->{help});
    die "$usage" unless ($args->{R1} and  ($args->{hisat2index} or $args->{bowtie2index}) and $args->{gtffile} );

    if ($args->{extraoptionsfile}) {
	die "ERROR: Could not locate extraoptionsfile $args->{extraoptionsfile}" if (! -e $args->{extraoptionsfile});
	$args->{extraoptionsfile} = abs_path($args->{extraoptionsfile}) 
    }

    ### Starttime ###
    $date = `date`;
    chomp $date;
    print "Starting at $date\n";
    &runtime;

    ### Standard parameters ###
    $args->{scratchfolder} = abs_path($args->{scratchfolder});
    die "R1 Fasta file $args->{R1} not found\n" if (!-e $args->{R1});
    $args->{R1} = abs_path($args->{R1});
    $args->{pe} = 0;
    if ($args->{R2}) {
        die "R2 Fasta file $args->{R2} not found\n" if (!-e $args->{R2});
        $args->{R2} = abs_path($args->{R2});
        $args->{pe} = 1;
    }
    $args->{logfile} = "$args->{scratchfolder}/log" unless ($args->{verbose});

    ### non-Standard parameters ###
    # add additional parameters processing here
    if ($args->{hisat2index}) {
	$args->{hisat2index} = hisat2indexcheck($args->{hisat2index});
    } else {
	$args->{bowtie2index} = bowtie2indexcheck($args->{bowtie2index});
    }
    die "ERROR: Could not locate gtffile $args->{gtffile}\n" if (! -e $args->{gtffile});
    $args->{gtffile} = abs_path($args->{gtffile});
    if ($args->{maskfile}) {
	die "ERROR: Could not locate maskfile $args->{maskfile}\n" if (! -e $args->{maskfile});
	$args->{maskfile} = abs_path($args->{maskfile});
    }

    ### Finish setup
    # check that required programs are installed
    requiredprograms(("fastqc", "\$TRIMMOMATIC", "hisat2", "cufflinks", "featureCounts", "\$PICARD"));

    diegracefully();
    readinextraoptions($args);
    scratchfoldersetup($args);

}
