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

Pipeline.pm - A Perl module providing functions for performing a variety of sequence analysis pipeline "stages"

=head1 DESCRIPTION

Put this module file in the same folder as a script using the module, then put these three lines in the script:

    use FindBin;                  # locate this script directory
    use lib "$FindBin::RealBin";  # look for modules in the script directory
    use Pipelinev3;               # load the pipeline module

=cut

######################################################

package Pipelinev4;
use Getopt::Long;
use Pod::Usage;
use Cwd 'abs_path';
use File::Temp qw( tempdir );
use File::Basename;
use feature "state";
use warnings;
use Exporter;
use Scalar::Util qw(looks_like_number);
use lib "$ENV{GOPHER_PIPELINES}/software/gbs-scripts"; # so we can load enzymes
use enzymes; # so we can assign padding sequences during samplesheet load

use IPC::Open3;

our @ISA= qw( Exporter );

# these functions CAN be exported.
# update this list with: grep "^sub" Pipelinev3.pm | cut -f2 -d' '
our @EXPORT_OK = qw( subsample subsampleplot allsubsample prinseq sizeselect fastqc qc trimmomatic trimmomaticplot kallisto kallistoplot salmon salmonplot hisat2 hisat2indexcheck hisat2plot tophat2 tophat2indexcheck tophat2plot bowtie2 bowtie2indexcheck bwa bwaindexcheck samtoolssortandindexbam samtoolsaddremoverg cleansam collectmultiplemetrics alignmentsummarymetrics alignmentsummarymetricsplot insertsize insertsizeplot wgsmetrics markduplicates hsmetrics hsmetricsplot hsmetricsprep samtoolsviewfilter samstat samtoolsstats featurecounts subreadplot cuffquant gbstrim allgbstrim allgunzip kraken krakendbcheck krakenplot bam2fastq freebayes freebayesplot vcfmetricsplot scratchfoldersetup writestats samplesheetandfoldersetup startreport getsinglesamplestats compilefolder fastqcplot expressiontableplot metrics allmetrics sphinx fixsamplesheet diegracefully readinextraoptions round10 round100 printlog run runtime formattime compressed INT_handler nodememory requiredprograms programversions title subtitle h1 h2 h3 h line fieldlist bullet space image link statprint $SCRATCHLOCAL $SCRATCHGLOBAL $GPRESOURCES);

# these functions are exported by default.
our @EXPORT = @EXPORT_OK;

#$SCRATCHLOCAL = "/scratch.local/"; # some scratches are full
$SCRATCHLOCAL = "/panfs/roc/scratch/";
$SCRATCHGLOBAL = "/panfs/roc/scratch/";
$GPRESOURCES = "$ENV{GOPHER_PIPELINES}/resources"; # gopher-pipelines resources directory

###########################################################################
###################### Single-sample pipeline stages ######################
###########################################################################

############################# Subsample #################################
# Requires: subsample.pl
# Args: subsample
# Stats: #rawreads #subsampledreads 
sub subsample {
    print "Running Subsample\n";

    my ($args, $stats) = @_;

    # determine sequence counts
    my $line = `wc -l $args->{R1}`;
    chomp $line;
    my ($count, $file) = split / /, $line;
    $stats->{subsample}{"#rawreads"} = $count / 4;
    $stats->{subsample}{"#subsampledreads"} = $count / 4;

    # skip subsampling if not called for
    if ((! $args->{subsample}) or $args->{subsample} == 0) {
        print "Skipping subsampling\n";
        return;
    }
    print "Subsampling $args->{subsample} reads\n";

    # subsample
    $newfile = $args->{R1};
    $newfile =~ s/\.fastq$/.subsample.fastq/;
    if ($args->{subsample} > $stats->{subsample}{"#rawreads"}) {
	`ln -s $args->{R1} $newfile`;
    } else {
	&run("subsample.pl $args->{R1} $args->{subsample} $count > $newfile", $args->{logfile});
    }
    $args->{R1} = $newfile;

    if ($args->{pe}) {
        $newfile = $args->{R2};
        $newfile =~ s/\.fastq$/.subsample.fastq/;
	if ($args->{subsample} > $stats->{subsample}{"#rawreads"}) {
	    `ln -s $args->{R2} $newfile`;
	} else {
	    &run("subsample.pl $args->{R2} $args->{subsample} $count > $newfile", $args->{logfile});
	}
        $args->{R2} = $newfile;
    }
    if ($args->{subsample} < $stats->{subsample}{"#rawreads"}) {
	$stats->{subsample}{"#subsampledreads"} = $args->{subsample};
    } else {
	$stats->{subsample}{"#subsampledreads"} = $stats->{subsample}{"#rawreads"};
    }

}

############################# All Subsample #################################
sub allsubsample {
    print "Running Subsample on all samples\n";

    my ($samples, $args, $stats) = @_;

    push @{$args->{stageorder}}, "subsample";

    # determine sequence counts
    my $myscratchfolder = "$args->{scratchfolder}/fastq-subsampled";
    mkdir $myscratchfolder;
    $gpout = "$myscratchfolder/gp.out";
    $progress = $args->{verbose} ? "--progress" : "";
    open GP, "| parallel $progress -j $args->{threads} > $gpout" or die "Cannot open pipe to gnu parallel\n";
#    open GP, ">$myscratchfolder/parallel.test" or die "Cannot open pipe to gnu parallel\n";

    # send commands to gnu parallel
    foreach $sample (keys %{$samples}) {
	print GP "wc -l $samples->{$sample}{R1}{fastq} >> $myscratchfolder/wc.out\n";
    }
    close GP;
    print "error code: $?\n" if ($?);
    open IFILE, "$myscratchfolder/wc.out" or die "Cannot open file $gpout: $!\n";
    $samplecount = 0;
    $sum = 0;
    while ($line = <IFILE>) {
	chomp $line;
	my ($count, $file) = split / /, $line;
	my ($fname, $path) = fileparse($file);
	my ($sample, $junk) = split /_R[12].fastq/, $fname;
	$stats->{$sample}{subsample}{"#rawreads"} = $count / 4;
	$stats->{$sample}{subsample}{"#subsampledreads"} = $count / 4;
	$sum += $count / 4;
	$samplecount++;
    }
    close IFILE;

    if ($args->{subsample} eq "auto") {
	$mean = int($sum / $samplecount);
	$max = $mean * 2;
	$min = $mean * .1; # need to do something with this
	print "Setting subsample to two times the mean = $max\n";
	$args->{subsample} = $max;
    }

    # skip subsampling if not called for
    if ($args->{subsample} == 0) {
	print "Skipping subsampling\n";
	return;
    }
    print "Subsampling $args->{subsample} reads from each sample\n";
    
    # run gnu parallel
    open GP, "| parallel $progress -j $args->{threads} > $gpout" or die "Cannot open pipe to gnu parallel\n";
#    open GP, ">parallel.test" or die "Cannot open pipe to gnu parallel\n";

    # send commands to gnu parallel
    foreach $sample (keys %{$samples}) {
	my ($fname, $path) = fileparse($samples->{$sample}{R1}{fastq});
	$newfile = "$myscratchfolder/$fname";
	# rich man's speedy subsampling using subsampler.pl
	$lines = $stats->{$sample}{subsample}{"#rawreads"} * 4; # four lines per read
	if ($args->{subsample} > $stats->{$sample}{subsample}{"#rawreads"}) {
	    `ln -s $samples->{$sample}{R1}{fastq} $newfile`;
	} else {
	    print GP "subsample.pl $samples->{$sample}{R1}{fastq} $args->{subsample} $lines > $newfile\n";
	}
	$samples->{$sample}{R1}{fastq} = $newfile;
	if ($samples->{$sample}{R2}{fastq}) {
	    ($fname, $path) = fileparse($samples->{$sample}{R2}{fastq});
	    $newfile = "$myscratchfolder/$fname";
	    # rich man's speedy subsampling using subsampler.pl
	    if ($args->{subsample} > $stats->{$sample}{subsample}{"#rawreads"}) {
		`ln -s $samples->{$sample}{R2}{fastq} $newfile`;
	    } else {
		print GP "subsample.pl $samples->{$sample}{R2}{fastq} $args->{subsample} $lines > $newfile\n";
	    }
	    $samples->{$sample}{R2}{fastq} = $newfile;	    
	}
	if ($args->{subsample} < $stats->{$sample}{subsample}{"#rawreads"}) {
	    $stats->{$sample}{subsample}{"#subsampledreads"} = $args->{subsample};
	} else {
	    $stats->{$sample}{subsample}{"#subsampledreads"} = $stats->{$sample}{subsample}{"#rawreads"};
	}
    }
    close GP;
    print "error code: $?\n" if ($?);
    $args->{fastqfolder} = "$myscratchfolder";

}

############################# Subsampleplot #################################
# Requires: ggplot2
sub subsampleplot {
    print "Running Subsampleplot\n";

    my ($samples, $args, $stats) = @_;

    # generate reads per sample plot
    $ofile = "subsampleplot.tmp";
    $rfile = "subsampleplot.r";
    print "Generating reads per sample plot\n" if ($args->{verbose});
    open OFILE, ">$ofile" or die "cannot open temporary file $ofile for writing: $!\n"; 
    # print out the data
    print OFILE "sample\tdata\tvalue\n";
    foreach $sample (sort keys %{$samples}) {
	$ss = $stats->{$sample}{subsample}{"#subsampledreads"} // 0;
	$raw = $stats->{$sample}{subsample}{"#rawreads"} - $ss;
	print OFILE "$sample\tSubsampledReads\t$ss\n";
	print OFILE "$sample\tRawReads\t$raw\n";
    }
    close OFILE;

    my $height = 480;
    my $width = 480;
    my $numsamples = keys %{$samples};
    if ($numsamples > 6) {
	$width = 480 + (20 * ($numsamples-6));
    }

    if ($args->{subsample}) {

    open RFILE, ">$rfile" or die "Cannot open $rfile\n";
    print RFILE qq(
  library(ggplot2);
  library(scales);
  datat <- read.table("$ofile", header=T, colClasses=c("sample"="factor"));
  png(filename="reads.v1.png", height = $height, width = $width);

  ggplot(datat, aes(x = sample, y = value, fill = data)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) + xlab("Sample") + ylab("Number of Reads") + theme(legend.title=element_blank()) + scale_y_continuous(labels = comma) + scale_fill_manual(values=c("#7F7F7F", "#00008B"));

  dev.off();

  png(filename="reads.v2.png", height = $height, width = $width);

  ggplot(datat, aes(reorder(sample, -value), value, fill=data)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) + xlab("Sample (sorted by Y-value)") + ylab("Number of Reads") + theme(legend.title=element_blank()) + scale_y_continuous(labels = comma) +scale_fill_manual(values=c("#7F7F7F", "#00008B"))

  dev.off();

  #eof
);

    close RFILE;
    } else {

    open RFILE, ">$rfile" or die "Cannot open $rfile\n";
    print RFILE qq(
  library(ggplot2);
  library(scales);
  datat <- read.table("$ofile", header=T, colClasses=c("sample"="factor"));
  png(filename="reads.v1.png", height = $height, width = $width);

  ggplot(datat, aes(x = sample, y = value, fill = data)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) + xlab("Sample") + ylab("Number of Reads") + theme(legend.position="none") + scale_y_continuous(labels = comma) + scale_fill_manual(values=c("#7F7F7F", "#00008B"))

  dev.off();

  png(filename="reads.v2.png", height = $height, width = $width);

  ggplot(datat, aes(reorder(sample, -value), value, fill=data)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) + xlab("Sample (sorted by Y-value)") + ylab("Number of Reads") + theme(legend.position="none") + scale_y_continuous(labels = comma) + scale_fill_manual(values=c("#7F7F7F", "#00008B"))

  dev.off();

  #eof
);

    close RFILE;
    }
    system("R --no-restore --no-save --no-readline < $rfile > $rfile.out");

    `mv *.png $args->{outputfolder}`;

    # add to report
    &h1("Reads per sample Plots");  
    &image("reads.v1.png");
    &space;
}


############################# Prinseq #################################
# Requires: prinseq
# Args:
# Stats: ??? 
sub prinseq {
    print "Running prinseq\n";

    my ($args, $stats) = @_;

    if ($args->{pe}) {
	&run("prinseq-lite.pl -fastq $args->{R1} -fastq2 $args->{R2} -out_good $args->{scratchfolder}/prinseq -out_bad $args->{scratchfolder}/prinseq.bad -log $args->{scratchfolder}/prinseq.log -custom_params \"A 70%;T 70%;G 70%;C 70%\"", $args->{logfile});
	# give the fastq files a better name
	my $newr1 = $args->{R1};
	$newr1 =~ s/\.fastq$/.prinseq.fastq/;
	my $newr2 = $args->{R2};
	$newr2 =~ s/\.fastq$/.prinseq.fastq/;
	`mv $args->{scratchfolder}/prinseq_1.fastq $newr1`;
	`mv $args->{scratchfolder}/prinseq_2.fastq $newr2`;
	$args->{R1} = $newr1;
	$args->{R2} = $newr2;
    } else {
	&run("prinseq-lite.pl -fastq $args->{R1} -out_good $args->{scratchfolder}/prinseq -out_bad $args->{scratchfolder}/prinseq.bad -log $args->{scratchfolder}/prinseq.log -custom_params \"A 70%;T 70%;G 70%;C 70%\"", $args->{logfile});
	`mv $args->{scratchfolder}/prinseq_1.fastq $args->{scratchfolder}/R1.prinseq.fastq`;
	# give the fastq file a better name
	my $newr1 = $args->{R1};
	$newr1 =~ s/\.fastq$/.prinseq.fastq/;
	`mv $args->{scratchfolder}/prinseq.fastq $newr1`;
	$args->{R1} = $newr1;
    }

    # save stats
    if ($args->{pe}) {
	$result = `grep "Good sequences (pairs):" $args->{scratchfolder}/prinseq.log`;
	$result =~ s/,//g;
	@result = split ' ', $result;
	my $goodsequences = $result[6];
	$stats->{prinseq}{goodsequences} = $goodsequences;
	
	$result = `grep "Input sequences (file 1):" $args->{scratchfolder}/prinseq.log`;
	$result =~ s/,//g;
	@result = split ' ', $result;
	my $inputsequences = $result[7];
	$stats->{prinseq}{inputsequences} = $inputsequences;
	$stats->{prinseq}{"%goodsequences"} = round100($goodsequences / $inputsequences * 100) if ($inputsequences > 0);
    } else {
	$result = `grep "Good sequences:" $args->{scratchfolder}/prinseq.log`;
	$result =~ s/,//g;
	@result = split ' ', $result;
	my $goodsequences = $result[5];
	$stats->{prinseq}{goodsequences} = $goodsequences;
	
	$result = `grep "Input sequences (file 1):" $args->{scratchfolder}/prinseq.log`;
	$result =~ s/,//g;
	@result = split ' ', $result;
	my $inputsequences = $result[5];
	$stats->{prinseq}{inputsequences} = $inputsequences;
	$stats->{prinseq}{"%goodsequences"} = round100($goodsequences / $inputsequences * 100) if ($inputsequences > 0);

    }
}

############################# Size select #################################
# Requires: pear
# Args:
# Stats: ??? 
sub sizeselect {
    print "Running sizeselect\n";

    my ($args, $stats) = @_;

    $min = $args->{sizeselectmin} // 0;
    $max = $args->{sizeselectmax} // 1000;

    ### run pear to stitch reads, run through stitched reads saving IDs of reads with good and bad lengths
    $count = 0;
    $goodcount = 0;
    if ($args->{pe}) {
	&run("pear -f $args->{R1} -r $args->{R2} -o $args->{scratchfolder}/pear ", $args->{logfile});
	$ifile = "$args->{scratchfolder}/pear.assembled.fastq";
	$ofile = "$args->{scratchfolder}/pear.ids.txt";
	open IFILE, $ifile or die "Cannot open $ifile: $!\n";
	open OFILE, $ofile or die "Cannot open $ofile: $!\n";
	while ($id = <IFILE>) {
	    $count++;
	    my $seq = <IFILE>;
	    my $plus = <IFILE>;
	    my $qual = <IFILE>;
	    $length = length($seq) - 1; # minus one for newline character
	    if (($length >= $min) and ($length <= $max)) {
		print OFILE "$id\n";
		$goodcount++;
	    }
	}
	close IFILE;
	close OFILE;

	### pull out good reads from original fastq file
	$newr1 = "$args->{R1}";
	$newr2 = "$args->{R2}";
	&run("extract_fastq.pl $ofile $args->{R1} > $newr1", $args->{logfile});
	&run("extract_fastq.pl $ofile $args->{R2} > $newr2", $args->{logfile});
	$args->{R1} = $newr1;
	$args->{R2} = $newr2;

	# save stats
	if ($args->{pe}) {
	    $stats->{sizeselect}{goodsequences} = $goodcount;
	    $stats->{sizeselect}{"%goodsequences"} = round100($goodcount / $count * 100) if ($count > 0);
	}
	
    } else {
	# how to handle single reads???
    }

}

######################## FastQC ##########################
# Requires: fastqc
# Args: 
# Stats: %gc lastq30baseR1 lastq20baseR1 meanreadqualityR1
sub fastqc {
    print "Running FastQC\n";

    my ($args, $stats) = @_;
    my $stage = $_[2] // "fastqc";

    my $myscratchfolder = $args->{scratchfolder} . "/$stage";
    mkdir $myscratchfolder;
    my $nogroup = ($args->{pacbio}) ? "" : "--nogroup";

    if ($args->{pe}) {
        &run("fastqc --threads $args->{threads} $nogroup -o $myscratchfolder $args->{R1} $args->{R2}", $args->{logfile});
    } else {
	&run("fastqc $nogroup -o $myscratchfolder $args->{R1}", $args->{logfile});
    }

    # Save FastQC stats
    # read in stats
    @reads = ("R1");
    push @reads, "R2" if ($args->{pe});
    foreach $read (@reads) {
#	($sample) = split /_/, $args->{$read};
	my ($filename, $path, $extension) = fileparse($args->{$read}, (".fastq"));
	$sample = $filename;
	$result = `unzip $myscratchfolder/${sample}_fastqc.zip */fastqc_data.txt -d $myscratchfolder`;
	print $result if ($args->{verbose});
	my $ifile = "$myscratchfolder/${sample}_fastqc/fastqc_data.txt";
	open IFILE, "$ifile" or die "Cannot open fastqc data file $ifile: $!\n";

	my $text;
	while ($line = <IFILE>) {
	    if ($line =~ /^>>Basic Statistics/) {
		while ($line = <IFILE>) {
		    last if ($line =~ /^>>END_MODULE/);
		    chomp $line;
		    if ($line =~ /^%GC/) {
			($text, $gc) = split /\t/, $line;
		    }
		    if ($line =~ /^Total Sequences/) {
			($text, $totalsequences) = split /\t/, $line;
		    }
		    if ($line =~ /^Sequence length/) {
			($text, $sequencelength) = split /\t/, $line;
		    }
		}
	    }		    
	    if ($line =~ /^>>Per base sequence quality/) {
		$header = <IFILE>;
		@mean = ();
		while ($line = <IFILE>) {
		    last if ($line =~ /^>>END_MODULE/);
		    chomp $line;
		    ($base, $mean) = split /\t/, $line;
		    ($base, $junk) = split /-/, $base if ($base =~ /-/); # for pacbio support
		    $mean[$base] = $mean;
		}
	    }
	    if ($line =~ /^>>Per sequence quality scores/) {
		$header = <IFILE>;
		$countsum = 0;
		$qualitysum = 0;
		while ($line = <IFILE>) {
		    last if ($line =~ /^>>END_MODULE/);
		    chomp $line;
		    ($quality, $count) = split /\t/, $line;
		    $qualitysum += $quality * $count;
		    $countsum += $count;
		}
	    }
	    if ($line =~ /^#Total Deduplicated Percentage/) {
		chomp $line;
		($junk, $deduppct) = split /\t/, $line;
	    }
	    if ($line =~ /^>>Adapter Content/) {
		my $oldline;
		while ($line = <IFILE>) {
		    if ($line =~ /^20\t/) {
			chomp $line;
			@line = split /\t/, $line;
			shift @line;
			$dimersum = 0;
			$dimersum += $_ for @line;
		    }
		    last if ($line =~ /^>>END_MODULE/);
		    $oldline = $line;
		}
		chomp $oldline;
		@line = split /\t/, $oldline;
		shift @line;
		$adaptersum = 0;
		$adaptersum += $_ for @line;
	    }
	}
	close IFILE;
	
	# calculate stats
	$q30 = 0;
	$q20 = 0;
	for $base (1..$#mean) {
	    $q30 = $base if ($mean[$base] >= 30);
	    $q20 = $base if ($mean[$base] >= 20);
	}
	$stats->{$stage}{"%gc"} = $gc;
	$stats->{$stage}{"totalsequences"} = $totalsequences;
	$stats->{$stage}{"sequencelength"} = $sequencelength;
	$stats->{$stage}{"lastq30base$read"} = $q30;
	$stats->{$stage}{"lastq20base$read"} = $q20;
	$stats->{$stage}{"%dimer"} = round100($dimersum);	
	$stats->{$stage}{"%adapter"} = round100($adaptersum);	
	$stats->{$stage}{"%deduplicated"} = round100($deduppct);
	$stats->{$stage}{"meanreadquality$read"} = round10($qualitysum / $countsum) if ($countsum > 0); # don't crash if dividing by zero
    }

}

############################# QC #################################
# Calculate some quality control statistics on a fastq file. Use of this 
# stage is deprecated, use the fastqc stage instead, which pulls many 
# of the same statistics from the fastqc output.
# Require: fastqQC.pl
# Args: 
# Stats: 
sub qc {
    print "Running QC\n";

    my ($args, $stats) = @_;

    my $ofile = "$args->{scratchfolder}/fastqQC.log";

    if ($args->{pe}) {
        $result = `fastqQC.pl $args->{R1} > $ofile.R1 & fastqQC.pl $args->{R2} > $ofile.R2; wait`;
        print $result;
    } else {
	$result = `fastqQC.pl $args->{R1} > $ofile.R1`;
	print $result;
    }

# Todo (low priority): save stats generated by fastqQC.pl into %stats
}

######################## Trimmomatic ##########################
# Run trimmomatic, supports extraoptions
# Require: java $TRIMMOMATIC
# Args: adapterfile (optional)
# Stats: inputreads surviving dropped %surviving %dropped
# Stats PE: inputreadpairs bothsurviving forwardonlysurvivng reverseonlysurviving dropped %bothsurviving %forwardonlysurvivng %reverseonlysurviving %dropped 
sub trimmomatic {
    print "Running trimmomatic\n";

    my ($args, $stats) = @_;

    # set trimmomatic parameters
    $xmx = "-Xmx2000M"; # java memory
    my $adapterfile = $args->{adapterfile} // "/panfs/roc/itascasoft/trimmomatic/0.33/adapters/all_illumina_adapters.fa";
    my $minscore = 16; # sliding window minimum q-score
    # set minimum length after trimming to be half of the original read length
    my $length = `sed '2q;d' $args->{R1} | wc -c`;
    $length--; # subtract newline
    my $minlength = int($length / 2);
    my $headcrop = $args->{headcrop} ? "HEADCROP:$args->{headcrop}" : "";

    my $trimmomaticoptions = $args->{extraoptions}{trimmomatic} // "ILLUMINACLIP:$adapterfile:2:30:10:2:true $headcrop LEADING:3 TRAILING:3 SLIDINGWINDOW:4:$minscore MINLEN:$minlength";

    $newr1 = $args->{R1};
    $newr1 =~ s/\.fastq$/.trim.fastq/;
    $newr1single = $newr1;
    $newr1single =~ s/\.fastq$/.singleton.fastq/;
    if ($args->{pe}) {
        $newr2 = $args->{R2};
        $newr2 =~ s/\.fastq$/.trim.fastq/;
        $newr2single = $newr2;
        $newr2single =~ s/\.fastq$/.singleton.fastq/;

        $result = &run("java $xmx -jar \$TRIMMOMATIC/trimmomatic.jar PE -phred33 -threads $args->{threads} $args->{R1} $args->{R2} $newr1 $newr1single $newr2 $newr2single $trimmomaticoptions", $args->{logfile});
    } else {
        $result = &run("java $xmx -jar \$TRIMMOMATIC/trimmomatic.jar SE -threads $args->{threads} $args->{R1} $newr1 $trimmomaticoptions", $args->{logfile});
    }

    # save stdout with all of the stats to a file
    my $ofile = "$args->{scratchfolder}/trimmomatic.log";
    open OFILE, ">$ofile" or print "Cannot open trimmomatic log $ofile: $!\n";
    print OFILE $result;
    close OFILE;

    # set trimmed fastq files as default fastq files
    $args->{R1} = $newr1;
    $args->{R2} = $newr2 if ($args->{pe});
}

### Trimmomaticplot ###
sub trimmomaticplot {
    my ($samples, $args, $stats) = @_;

    print "Running trimmomaticplot\n" if ($args->{verbose});
    open OFILE, ">trimmomatic-filelist.txt" or die "cannot open trimmomatic-filelist.txt: $!\n";
    foreach $sample (sort keys %{$samples}) {
            print OFILE "$args->{scratchfolder}/singlesamples/$sample/trimmomatic.log\t$sample\n";
    }
    close OFILE;
    @result = `trimmomaticplot.pl -f trimmomatic-filelist.txt`;

    # process the stdout, grabbing some stats
    $statzone = 0;
    foreach $result (@result) {
        $statzone = 0 if ($result =~ /^STATS End/);
        if ($statzone) {
            chomp $result;
            ($sample, $key, $stat) = split /\t/, $result;
            $stats->{$sample}{trimmomatic}{$key} = $stat;
        }
        $statzone = 1 if ($result =~ /^STATS Start/);
    }
    `mv *.png $args->{outputfolder}`;

    # add to report
    &h1("Trimmomatic Plots");  
    &image("trimmomatic.png");
    &image("trimmomatic-pct.png");
    &space;
}

######################## Kallisto ##########################
# Run RNA-Seq program Kallisto against a Kallisto index, supports extraoptions
# Require: kallisto samtools
# Args: kallistoindex
# Stats: 
sub kallisto {
    print "Running Kallisto\n";

    my ($args, $stats) = @_;
    my $myscratchfolder = $args->{scratchfolder} . "/kallisto";
    mkdir $myscratchfolder;

    my $stranded = $args->{stranded} ? "--rf-stranded" : "";
    my $kallistooptions = $args->{extraoptions}{kallisto} // "--bias --bootstrap-samples=100 $stranded"; # Set default settings here

    if ($args->{pe}) {
#	&run("kallisto quant $kallistooptions --threads $args->{threads} --index=$args->{kallistoindex} --pseduobam --output-dir=$myscratchfolder $args->{R1} $args->{R2} 2>>$args->{logfile} | samtools view -Sb - > $myscratchfolder/kallisto.bam", $args->{logfile}); # this generates a pseudobam
	$return = &run("kallisto quant $kallistooptions --threads $args->{threads} --index=$args->{kallistoindex} --output-dir=$myscratchfolder $args->{R1} $args->{R2}", $args->{logfile});
    } else {
	$return = &run("kallisto quant $kallistooptions --threads $args->{threads} --index=$args->{kallistoindex} --output-dir=$myscratchfolder --single $args->{R1}", $args->{logfile});
    }

    # write stdout to file so alignment metrics can be pulled from them later
    my $ofile = "$myscratchfolder/kallisto.log";
    open OFILE, ">$ofile" or print "Cannot open file $ofile for writing: $!\n";
    print OFILE $return;
    close OFILE;

}

### Kallistoplot ###
sub kallistoplot {
    my ($samples, $args, $stats) = @_;

    print "Running kallistoplot\n" if ($args->{verbose});
    open OFILE, ">kallisto-filelist.txt" or die "cannot open kallisto-filelist.txt: $!\n";
    foreach $sample (sort keys %{$samples}) {
            print OFILE "$args->{scratchfolder}/singlesamples/$sample/kallisto/kallisto.log\t$sample\n";
    }
    close OFILE;
    @result = `kallistoplot.pl -f kallisto-filelist.txt`;

    # process the stdout, grabbing some stats
    $statzone = 0;
    foreach $result (@result) {
        $statzone = 0 if ($result =~ /^STATS End/);
        if ($statzone) {
            chomp $result;
            ($sample, $key, $stat) = split /\t/, $result;
            $stats->{$sample}{kallisto}{$key} = $stat;
        }
        $statzone = 1 if ($result =~ /^STATS Start/);
    }
    `mv *.png $args->{outputfolder}`;

    # add to report
    &h1("Kallisto Plots");
    &image("kallisto.png");
    &space;
}

######################## Salmon ##########################
# Run RNA-Seq program Salmon against a Salmon index, supports extraoptions
# Require: salmon
# Args: salmonindex
# Stats: 
sub salmon {
    print "Running Salmon\n";

    my ($args, $stats) = @_;
    my $myscratchfolder = $args->{scratchfolder} . "/salmon";
    mkdir $myscratchfolder;

    my $stranded;
    if ($args->{pe}) {
	$stranded = $args->{stranded} ? "--libType ISR" : "--libType IU";
    } else {
	$stranded = $args->{stranded} ? "--libType SR" : "--libType U";
    }
    my $salmonoptions = $args->{extraoptions}{salmon} // "$stranded -k 31 --biasCorrect"; # set default settings here

    if ($args->{pe}) {
	&run("salmon quant $salmonoptions --index $args->{salmonindex} --threads $args->{threads} --output $myscratchfolder --mates1 $args->{R1} --mates2 $args->{R2}", $args->{logfile});
    } else {
	&run("salmon quant $salmonoptions --index $args->{salmonindex} --threads $args->{threads} --output $myscratchfolder --unmatedReads $args->{R1}", $args->{logfile});
    }

}

### Salmonplot ###
sub salmonplot {
    my ($samples, $args, $stats) = @_;

    print "Running salmonplot\n" if ($args->{verbose});
    open OFILE, ">salmon-filelist.txt" or die "cannot open salmon-filelist.txt: $!\n";
    foreach $sample (sort keys %{$samples}) {
            print OFILE "$args->{scratchfolder}/singlesamples/$sample/salmon/logs/salmon_quant.log\t$sample\n";
    }
    close OFILE;
    @result = `salmonplot.pl -f salmon-filelist.txt`;

    # process the stdout, grabbing some stats
    $statzone = 0;
    foreach $result (@result) {
        $statzone = 0 if ($result =~ /^STATS End/);
        if ($statzone) {
            chomp $result;
            ($sample, $key, $stat) = split /\t/, $result;
            $stats->{$sample}{salmon}{$key} = $stat;
        }
        $statzone = 1 if ($result =~ /^STATS Start/);
    }
    `mv *.png $args->{outputfolder}`;

    # add to report
    &h1("Salmon Plots");
    &image("salmon.png");
    &space;
}

######################## Hisat2 ##########################
# Align using hisat2
# Require: hisat2 samtools
# Args: hisat2index
# Stats: Many, see code for details, may just want to run picard alignment summmary
sub hisat2 {
    print "Running Hisat2\n";

    my ($args, $stats) = @_;

    my $bam = "$args->{scratchfolder}/hisat2.bam";
    my $stranded;

    if ($args->{pe}) {
	$stranded = $args->{stranded} ? "--rna-strandness RF" : "";
    } else {
	$stranded = $args->{stranded} ? "--rna-strandness R" : "";
    }
    my $hisat2options = $args->{extraoptions}{hisat2} // "--dta-cufflinks $stranded";

    # make tmp folder unique
    $srttmp = tempdir("sortXXXX", DIR => "$args->{scratchfolder}/", CLEANUP => 1);
    $srttmp .= "/tmp";

    # figure out how much ram to throw at samtools sort
    my $memory = nodememory();
    my $mem = int($memory / $args->{threads} * .8); # use memory proportional to the number of threads we're supposed to use, then reduce by 20% safety factor
    $mem = "${mem}G";

    my $result = "";
    if ($args->{pe}) {
#        $result = run("hisat2 $hisat2options --threads $args->{threads} -x $args->{hisat2index} -1 $args->{R1} -2 $args->{R2} 2>$args->{scratchfolder}/align_summary.txt | samtools sort -m $mem -T $srttmp -\@ $args->{threads} -O bam -o $bam", $args->{logfile});
        $result = run("hisat2 $hisat2options --threads $args->{threads} -x $args->{hisat2index} -1 $args->{R1} -2 $args->{R2} -S $bam", $args->{logfile});
    } else {
#        $result = run("hisat2 $hisat2options --threads $args->{threads} -x $args->{hisat2index} -U $args->{R1} 2>$args->{scratchfolder}/align_summary.txt | samtools sort -m $mem -T $srttmp -\@ $args->{threads} -O bam -o $bam", $args->{logfile});
        $result = run("hisat2 $hisat2options --threads $args->{threads} -x $args->{hisat2index} -U $args->{R1} -S $bam", $args->{logfile});
    }

    # send hisat output somewhere useful
    $ofile = "$args->{scratchfolder}/align_summary.txt";
    open OFILE, ">$ofile" or print "Cannot open file $ofile: $!\n";
    @results = split /\n/, $result;
    foreach $result (@results) {
	next if ($result =~ /^Warning/);
	print OFILE "$result\n";
    }
    close OFILE;

#    `echo $result > $args->{scratchfolder}/align_summary.txt`;

    die "Hisat2 failure\n" unless (-e $bam);

    # index bam
#    run("samtools index $bam", $args->{logfile});
    $args->{bam} = $bam;
    $stats->{files}{bam} = "$bam";
}

### verify that a hisat2 index is at the specified location, make it's path absolute
sub hisat2indexcheck {
    my ($index) = @_;
    my $extension = ".1.ht2";
    my $indexfile = $index . $extension;
    my $extension2 = ".1.ht2l";
    my $indexfile2 = $index . $extension2;
    if (-e $indexfile) {
	$indexfile = abs_path($indexfile);
    } elsif (-e $indexfile2) {
	$indexfile = abs_path($indexfile2);
    } else {
	die "Unable to find hisat2index $index";
    }
    ($name, $path, $suffix) = fileparse($indexfile, ($extension, $extension2));
    return "$path/$name";
}

### Hisat2plot ###
sub hisat2plot {
    my ($samples, $args, $stats) = @_;

    print "Running hisat2plot\n" if ($args->{verbose});
    open OFILE, ">hisat-filelist.txt" or die "cannot open hisat-filelist.txt: $!\n";
    foreach $sample (sort keys %{$samples}) {
#       if (-e "$scratchfolder/rnaseqQC-ss/$sample/tophat_out/spike.bam") {
#           print OFILE "$scratchfolder/rnaseqQC-ss/$sample/tophat_out/align_summary.txt\t$sample\trnaseqQC-ss/$sample/tophat_out/spike.sam\n";
#       } else {
            print OFILE "$args->{scratchfolder}/singlesamples/$sample/align_summary.txt\t$sample\n";
#       }
    }
    close OFILE;
    @result = `hisat2plot.pl -f hisat-filelist.txt`;

    # process the stdout, grabbing some stats
    $statzone = 0;
    foreach $result (@result) {
        $statzone = 0 if ($result =~ /^STATS End/);
        if ($statzone) {
            chomp $result;
            ($sample, $key, $stat) = split /\t/, $result;
            $stats->{$sample}{hisat2}{$key} = $stat;
        }
        $statzone = 1 if ($result =~ /^STATS Start/);
    }
    `mv *.png $args->{outputfolder}`;

    # add to report
    &h1("Hisat2 Plots");
    if ($args->{pe}) {
        &image("alignments-pairs.png");
        &image("alignments-pairs-pct.png");
    } else {
        &image("alignments.png");
        &image("alignments-pct.png");
    }
    &space;
}

######################## Tophat2 ##########################
# Align using tophat2
# Require: tophat2 samtools
# Args: tophat2index, gtffile
# Stats: Many, see code for details, may just want to run picard alignment summmary
sub tophat2 {
    print "Running Tophat2\n";

    my ($args, $stats) = @_;

    my $stranded = ""; # stranded alignment not supported yet
    if ($args->{pe}) {
	$stranded = $args->{stranded} ? "--library-type fr-firststrand" : "";
    } else {
	$stranded = $args->{stranded} ? "--library-type fr-firststrand" : "";
    }

    my $tophat2options = $args->{extraoptions}{tophat2} // "--no-coverage-search $stranded --GTF $args->{gtffile}";

    if ($args->{pe}) {
        # get insert size and standard deviation
        $result = run("insertsize.pl -p $args->{threads} -m .1 $args->{bowtie2index} $args->{R1} $args->{R2}", $args->{logfile});
        chomp $result;
        my ($mean, $stddev) = split /\t/, $result;
        $mean = int($mean);
        $stddev = int($stddev);
        `echo "Insert size mean: $mean" >> $args->{logfile}`;
        `echo "Insert size standard deviation: $stddev" >> $args->{logfile}`;
        my $readlength = `fastqreadlength.pl $args->{R1}`;
        chomp $readlength;
        my $innerdist = int($mean - ($readlength * 2));
        run("tophat $tophat2options -p $args->{threads} -r $innerdist --mate-std-dev $stddev -o $args->{scratchfolder}/tophat_out $args->{bowtie2index} $args->{R1} $args->{R2}", $args->{logfile});
    } else {
        run("tophat $tophat2options -p $args->{threads} -o $args->{scratchfolder}/tophat_out $args->{bowtie2index} $args->{R1}", $args->{logfile});
    }

    if (0) { # no longer need to sort here
    # make tmp folder unique
    $srttmp = tempdir("sortXXXX", DIR => "$args->{scratchlocal}/", CLEANUP => 1);
    $srttmp .= "/tmp";

    # figure out how much ram to throw at samtools sort
    my $memory = nodememory();
    my $mem = int($memory / $args->{threads} * .8); # use memory proportional to the number of threads we're supposed to use, then reduce by 20% safety factor
    $mem = "${mem}G";

    ### sort and index
    run("mv $args->{scratchfolder}/tophat_out/accepted_hits.bam $args->{scratchfolder}/tophat_out/accepted_hits.unsorted.bam", $args->{logfile});
    run("samtools sort -O bam -m $mem -T $srttmp $args->{scratchfolder}/tophat_out/accepted_hits.unsorted.bam > $args->{scratchfolder}/tophat_out/accepted_hits.bam", $args->{logfile});
    run("samtools index $args->{scratchfolder}/tophat_out/accepted_hits.bam", $args->{logfile});

    run("mv $args->{scratchfolder}/tophat_out/unmapped.bam $args->{scratchfolder}/tophat_out/unmapped.unsorted.bam", $args->{logfile});
    run("samtools sort -O bam -T $srttmp $args->{scratchfolder}/tophat_out/unmapped.unsorted.bam > $args->{scratchfolder}/tophat_out/unmapped.bam", $args->{logfile});
    run("samtools index $args->{scratchfolder}/tophat_out/unmapped.bam", $args->{logfile});
    }

    my $bam = "$args->{scratchfolder}/tophat_out/accepted_hits.bam";
    die "Tophat2 failure\n" unless (-e $bam);

    # index bam
    run("samtools index $bam", $args->{logfile});
    $args->{bam} = $bam;
    $stats->{files}{bam} = "$bam";
}

### verify that a tophat2 index is at the specified location, make it's path absolute
sub tophat2indexcheck {
    my ($index) = @_;
    my $extension = ".1.bt2";
    my $indexfile = $index . $extension;
    my $extension2 = ".1.bt2l";
    my $indexfile2 = $index . $extension2;
    if (-e $indexfile) {
	$indexfile = abs_path($indexfile);
    } elsif (-e $indexfile2) {
	$indexfile = abs_path($indexfile2);
    } else {
	die "Unable to find tophat2index $index";
    }
    ($name, $path, $suffix) = fileparse($indexfile, ($extension, $extension2));
    return "$path/$name";
}

### Tophat2plot ###
sub tophat2plot {
    my ($samples, $args, $stats) = @_;

    print "Running tophat2plot\n" if ($args->{verbose});
    open OFILE, ">tophat-filelist.txt" or die "cannot open tophat-filelist.txt: $!\n";
    foreach $sample (sort keys %{$samples}) {
#       if (-e "$scratchfolder/rnaseqQC-ss/$sample/tophat_out/spike.bam") {
#           print OFILE "$scratchfolder/rnaseqQC-ss/$sample/tophat_out/align_summary.txt\t$sample\trnaseqQC-ss/$sample/tophat_out/spike.sam\n";
#       } else {
            print OFILE "$args->{scratchfolder}/singlesamples/$sample/tophat_out/align_summary.txt\t$sample\n";
#       }
    }
    close OFILE;
    @result = `tophatplot.pl -f tophat-filelist.txt`;

    # process the stdout, grabbing some stats
    $statzone = 0;
    foreach $result (@result) {
        $statzone = 0 if ($result =~ /^STATS End/);
        if ($statzone) {
            chomp $result;
            ($sample, $key, $stat) = split /\t/, $result;
            $stats->{$sample}{tophat2}{$key} = $stat;
        }
        $statzone = 1 if ($result =~ /^STATS Start/);
    }
    `mv *.png $args->{outputfolder}`;

    # add to report
    &h1("Tophat2 Plots");
    if ($args->{pe}) {
        &image("alignments-pairs.png");
        &image("alignments-pairs-pct.png");
    } else {
        &image("alignments.png");
        &image("alignments-pct.png");
    }
    &space;
}

######################## Bowtie2 ##########################
# Align using bowtie2, pipe output through samtools sort to produce an indexed bam, supports extraoptions
# Require: bowtie2 samtools
# Args: bowtie2index
# Stats: Many, see code for details, may just want to run picard alignment summmary
sub bowtie2 {
    print "Running Bowtie2\n";

    my ($args, $stats) = @_;

    my $bam = "$args->{scratchfolder}/bowtie2.bam";
    my $bowtieoptions = $args->{extraoptions}{bowtie2} // "";

    # make tmp folder unique
    my $srttmp = tempdir("sortXXXX", DIR => "$args->{scratchfolder}/", CLEANUP => 1);
    $srttmp .= "/tmp";

    # figure out how much ram to throw at samtools sort
    my $memory = nodememory();
    $mem = int($memory / $args->{threads} * .8); # use memory proportional to the number of threads we're supposed to use, then reduce by 20% safety factor
    $mem = "${mem}G";

    $readgroup = ($args->{readgroup}) ? "--rg-id $args->{readgroup} --rg SM:$args->{readgroup}" : "";

    ### Run bowtie2 and samtools sort
    # Build commands
    my $bowtiecmd;
    if ($args->{pe}) {
	$bowtiecmd = "bowtie2 $bowtieoptions $readgroup -p $args->{threads} -x $args->{bowtie2index} -1 $args->{R1} -2 $args->{R2}";
    } else {
	$bowtiecmd = "bowtie2 $bowtieoptions $readgroup -p $args->{threads} -x $args->{bowtie2index} -U $args->{R1}";
    }
    my $samtoolscmd = "samtools sort -m $mem -T $srttmp -\@ $args->{threads} -O bam -o $bam >> $args->{logfile}";

    # Be a tricksy hobbit and collect stderr and stdout from bowtie2. 
    # Forward the stdout (sam data) to samtools, and forward the stderr
    # to the log and save a copy of stderr so we can parse out the alignment
    # metrics.
    open SAMTOOLS, "| $samtoolscmd" or die "Cannot start samtools: $!\n";
    open LOG, ">>$args->{logfile}" or die "Cannot open log file $args->{logfile}: $!\n";
    my $pid = open3(\*WRITER, \*READER, \*ERROR, $bowtiecmd);
    $result = "";
    # send stdout to samtools 
    while( my $output = <READER> ) {
	print SAMTOOLS $output; 
    }
    # send stderr to logfile and local variable
    while( my $errout = <ERROR> ) {
	print LOG $errout;
	$result .= $errout;
    }
    waitpid( $pid, 0 ) or die "$!\n";
    my $retval =  $?;
    close SAMTOOLS;
    close LOG;

    die "Bowtie2 failure\n" unless (-e $bam);

    # index bam
    &run("samtools index $bam");
    $args->{bam} = "$bam";
    $stats->{files}{bam} = "$bam";

    # Parse bowtie alignment stats
    @rawlines = split /\n/, $result;
    # filter out Warning lines
    foreach $line (@rawlines) {
	next if ($line =~ /^Warning/);
	push @lines, $line;
    }
    
    if ($args->{pe}) {
	$lines[0] =~ /^(\d+)/;
	$stats->{bowtie2}{reads} = $1;
	$lines[2] =~ /^\s+(\d+)/;
	$stats->{bowtie2}{concordantzero} = $1;
	$lines[3] =~ /^\s+(\d+)/;
	$stats->{bowtie2}{concordantsingle} = $1;
	$lines[4] =~ /^\s+(\d+)/;
	$stats->{bowtie2}{concordantmulti} = $1;
	$lines[7] =~ /^\s+(\d+)/;
	$stats->{bowtie2}{discordantsingle} = $1;
	
    } else {
	$lines[0] =~ /^(\d+)/;
	$stats->{bowtie2}{reads} = $1;
	$lines[2] =~ /^\s+(\d+)/;
	$stats->{bowtie2}{alignedzero} = $1;
	$lines[3] =~ /^\s+(\d+)/;
	$stats->{bowtie2}{alignedsingle} = $1;
	$lines[4] =~ /^\s+(\d+)/;
	$stats->{bowtie2}{alignedmulti} = $1;
    }

}

# verify that a bowtie2 index is at the specified location, make it's path absolute
sub bowtie2indexcheck {
    my ($index) = @_;
    my $extension = ".1.bt2";
    my $indexfile = $index . $extension;
    my $extension2 = ".1.bt2l";
    my $indexfile2 = $index . $extension2;
    if (-e $indexfile) {
	$indexfile = abs_path($indexfile);
    } elsif (-e $indexfile2) {
	$indexfile = abs_path($indexfile2);
    } else {
	die "Unable to find bowtie2index $index";
    }
    ($name, $path, $suffix) = fileparse($indexfile, ($extension, $extension2));
    return "$path/$name";
}

######################## BWA ##########################
# Align using bwa, pipe output through samtools sort to produce an indexed bam
# Require: bwa samtools
# Args: bwaindex
# Stats: none, run picard alignment summary statistics instead
sub bwa {
    print "Running BWA\n";

    my ($args, $stats) = @_;

    my $bam = "$args->{scratchfolder}/bwa.bam";

    # make tmp folder unique 
    my $srttmp = tempdir("sortXXXX", DIR => "$args->{scratchfolder}/", CLEANUP => 1);
    $srttmp .= "/tmp";

    # figure out how much ram to throw at samtools sort
    my $memory = nodememory();
    my $mem = int($memory / $args->{threads} * .8); # use memory proportional to the number of threads we're supposed to use, then reduce by 20% safety factor
    $mem = "${mem}G";
#    print STDERR "samtools sort mem: $mem\n";
    $mem = $args->{bwamem} // $mem; 

    $readgroup = ($args->{readgroup}) ? "-R '\@RG\\tID:$args->{readgroup}\\tSM:$args->{readgroup}'" : "";

    my $logging = ($args->{verbose}) ? "" : "2>>$args->{logfile}";
    if ($args->{pe}) {
        &run("bwa mem -v 3 -t $args->{threads} $readgroup $args->{bwaindex} $args->{R1} $args->{R2} $logging | samtools sort -m $mem -T $srttmp -\@ $args->{threads} -O bam -o $bam 1>&2", $args->{logfile});
    } else {
	&run("bwa mem -t $args->{threads} $readgroup $args->{bwaindex} $args->{R1} $logging | samtools sort -m $mem -T $srttmp -\@ $args->{threads} -O bam -o $bam 1>&2", $args->{logfile});
    }
    die "BWA failure\n" unless (-e $bam);

    # index bam
    run("samtools index $bam");
    $args->{bam} = "$bam";
    $stats->{files}{bam} = "$bam";
}

# verify that a bwa index is at the specified location, make it's path absolute
sub bwaindexcheck {
    my ($index) = @_;
    my $extension = ".bwt";
    my $indexfile = $index . $extension;
    if (! -e $indexfile) {
	die "Unable to find bwaindex $index";
    }
    $indexfile = abs_path($indexfile);
    ($name, $path, $suffix) = fileparse($indexfile, ($extension));
    return "$path/$name";
}

########################### samtools sort and index #################
# sort and index a bam file
# Require: samtools
# Args: bam
# Stats: 
sub samtoolssortandindexbam {

    my ($args, $stats) = @_;
    print "Running Samtools sort and index\n";

    ### sort and index
    # make tmp folder unique
    my $srttmp = tempdir("sortXXXX", DIR => "$args->{scratchfolder}/", CLEANUP => 1);

    # figure out how much ram to throw at samtools sort
    my $memory = nodememory();
    my $mem = int($memory / $args->{threads} * .8); # use memory proportional to the number of threads we're supposed to use, then reduce by 20% safety factor
    $mem = "${mem}G";

    my $newbam = $args->{bam};
    $newbam =~ s/\.bam$/.sort.bam/;

    &run("samtools sort -O bam -m $mem -T ${srttmp}/ $args->{bam} -o $newbam", $args->{logfile});
    &run("samtools index $newbam", $args->{logfile});

    $args->{bam} = $newbam;
    die "Sort failure\n" if (! -e $args->{bam});
}

########################### Remove aligned #################
# Generate fastq file with aligned reads removed (for removing contaminants)
# Require: gopher-biotools, samtools
# Args: fastq, bam
# Stats: 
sub removealigned {

    my ($args, $stats) = @_;
    print "Running remove aligned\n";

    $idfile = "$args->{scratchfolder}/removealignedids.txt";

    # save aligned read ids to file
    $cleanedR1 = $args->{R1};
    $cleanedR1 =~ "$args->{scratchfolder}/R1.removealigned.fastq";
    &run("samtools view $args->{bam} | cut -f1 | sort | uniq > $idfile", $args->{logfile});
    # generate new fastq file with reads removed
    $result = &run("extract_fastq.pl -discard $idfile $args->{R1} > $cleanedR1", $args->{logfile});
    print $result if ($args->{verbose});
    $args->{R1} = $cleanedR1;

    if ($args->{R2}) {
	$cleanedR2 = $args->{R2};
	$cleanedR2 =~ "$args->{scratchfolder}/R2.removealigned.fastq";
	# generate new fastq file with reads removed
	$result = &run("extract_fastq.pl -discard $idfile $args->{R1} > $cleanedR2", $args->{logfile});
	print $result if ($args->{verbose});
	$args->{R2} = $cleanedR2;
    }

}


########################### coverage blocks summary #################
# Calculate blocks of coverage across the genome
# Require: bedtools, gbs-scripts
# Args: bam
# Stats: 
sub coverageblocks {

    # determine the number of loci with a minimum length of X bp and minimum depth of Y

    my ($args, $stats) = @_;
    print "Running coverage summary\n";

    &run("coverage-summary.pl --bamfile $args->{bam} --referencefai $args->{referencefasta}.fai", $args->{logfile});

}


########################### samtools add read group #################
# Add a read group to a bam file
# Require: samtools
# Args: bam
# Stats: 
sub samtoolsaddreplacerg {

    my ($args, $stats) = @_;
    print "Running Samtools add read group\n";

    ### add read group
    $tmpbam = $args->{bam} . ".tmp";
    `mv $args->{bam} $tmpbam`;
    &run("samtools addreplacerg -o $args->{bam} $tmpbam", $args->{logfile});

}

###################### Picard CleanSam ######################
# Requires: java, $PICARD
# Args: 
# Stats: 
sub cleansam {
    print "Running Picard CleanSam\n";

    my ($args, $stats) = @_;

    my $newbam = $args->{bam};
    $newbam =~ s/\.bam$/.clean.bam/;

    &run("java -Xmx4g -jar \$PICARD/picard.jar CleanSam INPUT=$args->{bam} OUTPUT=$newbam", $args->{logfile});

    $args->{bam} = $newbam;

}

###################### Picard CollectMultipleMetrics ######################
# Requires: java, $PICARD
# Args: referencefasta
# Stats: 
sub collectmultiplemetrics {
    print "Running Picard CollectMultipleMetrics\n";

    my ($args, $stats) = @_;

    if ($args->{pe}) {
	$defaultmetrics = "PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics";
    } else {
	$defaultmetrics = "PROGRAM=CollectAlignmentSummaryMetrics";
    }
    my $picardmodules = $_[2] // "PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics";

    &run("java -Xmx4g -jar \$PICARD/picard.jar CollectMultipleMetrics PROGRAM=null $picardmodules INPUT=$args->{bam} OUTPUT=$args->{scratchfolder}/picard REFERENCE_SEQUENCE=$args->{referencefasta}", $args->{logfile});

}

###################### Picard AlignmentSummaryMetrics ######################
# Requires: java, $PICARD
# Args: referencefasta
sub alignmentsummarymetrics {
    print "Running Picard AlignmentSummaryMetrics\n";

    my ($args, $stats) = @_;

    # load java module???
    &run("java -Xmx4g -jar \$PICARD/picard.jar CollectAlignmentSummaryMetrics INPUT=$args->{bam} OUTPUT=$args->{scratchfolder}/picard.alignment_summary_metrics REFERENCE_SEQUENCE=$args->{referencefasta}", $args->{logfile});

}

### picardalignmentsummarymetricsplot ###
sub alignmentsummarymetricsplot {
    my ($samples, $args, $stats) = @_;

    print "Gennerating picardalignmentsummarymetrics Plot\n" if ($args->{verbose});
    open OFILE, ">picardalignmentsummarymetrics-filelist.txt" or die "cannot open picardallignmentsummarymetrics-filelist.txt: $!\n";
    foreach $sample (sort keys %{$samples}) {
	print OFILE "$args->{scratchfolder}/singlesamples/$sample/picard.alignment_summary_metrics\t$sample\n";
    }
    close OFILE;
    @result = `picardalignmentsummarymetricsplot.pl -f picardalignmentsummarymetrics-filelist.txt`;
    # save stats from stdout
    $statzone = 0;
    foreach $result (@result) {
#        print $result if ($args->{verbose});
        $statzone = 0 if ($result =~ /^STATS End/);
        if ($statzone) {
            chomp $result;
            ($sample, $key, $stat) = split /\t/, $result;
            $stats->{$sample}{alignmentsummarymetrics}{$key} = $stat;
        }
        $statzone = 1 if ($result =~ /^STATS Start/);
    }
    `mv *.png $args->{outputfolder}`;

    # add to report
    &h1("Picard AlignmentSummaryMetrics Plots");
    if ($args->{pe}) {
	&image("alignments-pairs.v1.png");
	&image("alignments-pairs-pct.png");
    } else {
	&image("alignments.v1.png");
	&image("alignments-pct.png");
    }
    &space;

}

######################## Picard Insertsize Metrics ##########################
# Requires: java, $PICARD
sub insertsize {
    print "Running Picard insert size metrics\n";

    my ($args, $stats) = @_;

    unless ($args->{pe}) {
        print "Skipping Picard insert size metrics - not a paired-end dataset";
    }

    # load java module???
    &run("java -Xmx4g -jar \$PICARD/picard.jar CollectInsertSizeMetrics INPUT=$args->{bam} OUTPUT=$args->{scratchfolder}/picard.insert_size_metrics HISTOGRAM_FILE=$args->{scratchfolder}/picard.insert_size_histogram.pdf", $args->{logfile});
}

### insertsizeplot ###
sub insertsizeplot {
    my ($samples, $args, $stats, $plottype) = @_;

    if ($args->{pe}) {
	print "Generating insertsize plot\n" if ($args->{verbose});
	open OFILE, ">insertmetrics-filelist.txt" or die "cannot open insertmetrics-filelist.txt: $!\n";
	foreach $sample (sort keys %{$samples}) {
	    print OFILE "$args->{scratchfolder}/singlesamples/$sample/picard.insert_size_metrics\t$sample\n";
	}
	close OFILE;
	@result = `insertplot.pl -f insertmetrics-filelist.txt`;
	`mv insertplot.png $args->{outputfolder}`;
	`mv insertboxplot.png $args->{outputfolder}`;
	# save stats from stdout
	$statzone = 0;
	foreach $result (@result) {
	    $statzone = 0 if ($result =~ /^STATS End/);
	    if ($statzone) {
		chomp $result;
		($sample, $key, $stat) = split /\t/, $result;
		$stats->{$sample}{insertmetrics}{$key} = $stat;
	    }
	    $statzone = 1 if ($result =~ /^STATS Start/);
	}

	# add to report
	&h1("Picard InsertMetrics Plot");
	if ($plottype && $plottype eq "violin") {
	    &image("insertplot.png"); # don't show the violin plot any more
	} else {
	    &image("insertboxplot.png");
	}
	&space();
    }
}


###################### Picard CollectWgsMetrics ######################
# Requires: java, $PICARD
# Args: referencefasta
sub wgsmetrics {
    print "Running Picard CollectWgsMetrics\n";

    my ($args, $stats) = @_;

    # load java module???
    &run("java -Xmx4g -jar \$PICARD/picard.jar CollectWgsMetrics INPUT=$args->{bam} COVERAGE_CAP=null OUTPUT=$args->{scratchfolder}/picard.wgs_metrics REFERENCE_SEQUENCE=$args->{referencefasta}", $args->{logfile});

}

######################## Picard MarkDuplicates ##########################
# Requires: java, $PICARD
sub markduplicates {
    print "Running Picard MarkDuplicates\n";

    my ($args, $stats) = @_;

    my $newbam = $args->{bam};
    $newbam =~ s/\.bam$/.dedup.bam/;

    # load java module???
    &run("java -Xmx4g -jar \$PICARD/picard.jar MarkDuplicates INPUT=$args->{bam} OUTPUT=$newbam CREATE_INDEX=true METRICS_FILE=removeduplicates.txt", $args->{logfile});
    die "Mark Duplicates failure\n" unless (-e $newbam);
    $args->{bam} = $newbam;

    # grab some stats
    $ifile = "$args->{scratchfolder}/removeduplicates.txt";
    open IFILE, "$ifile" or die "Cannot open picard log file $ifile: $!\n";
    while ($line = <IFILE>) {
	if ($line =~ /## METRICS CLASS/) {
	    $line = <IFILE>;
	    chomp $line;
	    @keys = split /\t/, $line;
	    $line = <IFILE>;
	    chomp $line;
	    @values = split /\t/, $line;
	    for $i (0..$#keys) {
		$values[$i] = "" unless (defined($values[$i])); # fill in for empty values
		$data{$keys[$i]} = $values[$i];
	    }
	    last;
	}
    }

    $stats->{removeduplicates}{"\%duplication"} = $data{PERCENT_DUPLICATION};

}

######################## Picard HsMetrics ##########################
sub hsmetrics {
    print "Running Picard HsMetrics\n";

    my ($args, $stats) = @_;

    # load java module???
    $result = run("java -Xmx2g -jar \$PICARD/picard.jar CalculateHsMetrics INPUT=$args->{bam} OUTPUT=$args->{scratchfolder}/hsmetrics.txt BAIT_INTERVALS=$args->{baitfile} TARGET_INTERVALS=$args->{targetfile}", $args->{logfile});
# this version uses the reference genome to calculate more stats
#    $result = run("java -Xmx2g -jar \$PICARD/picard.jar CalculateHsMetrics INPUT=$args->{bam} OUTPUT=$args->{scratchfolder}/hsmetrics.txt BAIT_INTERVALS=$args->{scratchfolder}/baits.txt TARGET_INTERVALS=$args->{scratchfolder}/targets.txt REFERENCE_SEQUENCE=$args->{referencefasta} PER_TARGET_COVERAGE=$args->{scratchfolder}/pertargetcoverage.txt", $args->{logfile});

    print $result;
}

### HsMetrics Plot ###
sub hsmetricsplot {
    my ($samples, $args, $stats) = @_;

    print "Running hsmetricsplot\n" if ($args->{verbose});
    open OFILE, ">hsmetrics-filelist.txt" or die "cannot open hsmetrics-filelist.txt: $!\n";
    foreach $sample (sort keys %{$samples}) {
	print OFILE "$args->{scratchfolder}/singlesamples/$sample/hsmetrics.txt\t$sample\n";
    }
    close OFILE;
    @result = `picardhsmetricsplot.pl -f hsmetrics-filelist.txt`;

    # process the stdout, grabbing some stats
    $statzone = 0;
    foreach $result (@result) {
        $statzone = 0 if ($result =~ /^STATS End/);
        if ($statzone) {
            chomp $result;
            ($sample, $key, $stat) = split /\t/, $result;
            $stats->{$sample}{hsmetrics}{$key} = $stat;
        }
        $statzone = 1 if ($result =~ /^STATS Start/);
    }
    `mv *.png $args->{outputfolder}`;

    my $firstsample = (keys %{$stats})[0];

    # add to report
    &h1("Picard HsMetrics Plots");
    &fieldlist("Bait territory", "$stats->{$firstsample}{hsmetrics}{BAIT_TERRITORY}bp");
    &fieldlist("Target territory", "$stats->{$firstsample}{hsmetrics}{TARGET_TERRITORY}bp");
    &fieldlist("Genome size", "$stats->{$firstsample}{hsmetrics}{GENOME_SIZE}bp");
    &space;
    &image("picard-bait.png");
    &image("picard-coverage.png");
    &space;
}

### HsMetrics Prep ###
sub hsmetricsprep {
    my ($args) = @_;

    print "Running hsMetrics prep\n" if ($args->{verbose});

    # look for dict file in a couple obvious locations
    if (! $args->{referencedict}) {
	if ($args->{referencefasta}) {
	    my $tmpdict1 = $args->{referencefasta} . ".dict";
	    my $tmpdict2 = $args->{referencefasta};
	    if ($tmpdict2 =~ s/\.fa$/.dict/) {

	    } elsif ($tmpdict2 =~ s/\.fasta$/.dict/) {

	    } else {
		$tmpdict2 = "";
	    }
	    if (-e $tmpdict1) {
		$args->{referencedict} = $tmpdict1;
	    } elsif (-e $tmpdict2) {
		$args->{referencedict} = $tmpdict2;
	    }
	}
    }

    # create reference dict file if it doesn't already exist
    if (! $args->{referencedict}) {
	$args->{referencedict} = "$args->{scratchfolder}/reference.dict";
	run("java -Xmx2g -jar \$PICARD/picard.jar CreateSequenceDictionary REFERENCE=$args->{referencefasta} OUTPUT=$args->{referencedict}", $args->{logfile});
    }

    # convert BED to interval file
    $args->{baitfile} = "$args->{scratchfolder}/baits.interval_list";
    $args->{targetfile} = "$args->{scratchfolder}/targets.interval_list";

    run("java -Xmx2g -jar \$PICARD/picard.jar BedToIntervalList INPUT=$args->{baitbed} OUTPUT=$args->{baitfile} SEQUENCE_DICTIONARY=$args->{referencedict}", $args->{logfile});
    run("java -Xmx2g -jar \$PICARD/picard.jar BedToIntervalList INPUT=$args->{targetbed} OUTPUT=$args->{targetfile} SEQUENCE_DICTIONARY=$args->{referencedict}", $args->{logfile});

}

######################## Samtools View Filter ##########################
# Requires: samtools
sub samtoolsviewfilter {
    print "Running Samtools View filter\n";

    my ($args, $stats) = @_;

    my $newbam = $args->{bam};
    $newbam =~ s/\.bam$/.filter.bam/;

# include flags: read is paired: 1; read mapped in proper pair: 2 (total: 3)
# exclude flags: read is unmapped: 4; mate is unmapped: 8; read fails vendor QC: 512; read is PCR duplicate: 1024 (total: 1548)
    &run("samtools view -b -F 1548 -o $newbam $args->{bam}", $args->{logfile});
    die "Samtools view filter failure\n" unless (-e $newbam);
    &run("samtools index $newbam");
    $args->{bam} = $newbam;

}

######################## Samstat from 2013-07-08 ##########################
sub samstat {
    print "Running SAMStat\n";

    my ($args, $stats) = @_;

    run("samstat $args->{bam}", $args->{logfile});
    # give output file a better name
    `mv $args->{bam}.html $args->{scratchfolder}/samstat.html`;

}

######################## Samtools stats ##########################
sub samtoolsstats {
    print "Running Samtools stats plot\n";

    my ($args, $stats) = @_;

    my $ref = ($args->{referencefasta}) ? "--ref-seq $args->{referencefasta}" : "";
    my $statsfile = "$args->{scratchfolder}/samtools.stats";
    run("samtools stats $ref $args->{bam} > $statsfile", $args->{logfile});
    run("plot-bamstats -p $args->{scratchfolder}/samtoolsstats/ $statsfile", $args->{logfile});

    # get the peak indel number and %
    @lines = `grep ^IC $statsfile | cut -f 2-`;
    $maxsum = 0;
    foreach $line (@lines) {
	chomp $line;
	@line = split /\t/, $line;
	my $cycle = shift @line;
	$sum = $line[0] + $line[1] + $line[2] + $line[3];
	$maxsum = $sum if ($sum > $maxsum);	
    }

    $stats->{samtoolsstats}{peakindel} = $maxsum;
#    $stats->{samtoolsstats}{"%peakindel"} = $maxsum / $stats->{fastqc}{};

    # get the peak mismatch number and %
    @lines = `grep ^MPC $statsfile | cut -f 2-`;
    $maxsum = 0;
    foreach $line (@lines) {
	chomp $line;
	@line = split /\t/, $line;
	$sum = 0;
	for $index (31..$#line) {
	    $sum += $line[$index];
	}
	$maxsum = $sum if ($sum > $maxsum);
    }

    $stats->{samtoolsstats}{peakmismatch} = $maxsum;

}


######################## Subread featurecounts ##########################
sub featurecounts {
    print "Estimating abundances with Subread\n";

    my ($args, $stats) = @_;

    my $featurecountsoptions = $args->{extraoptions}{featurecounts} // "";

    run("featureCounts -T $args->{threads} $featurecountsoptions -a $args->{gtffile} -o $args->{scratchfolder}/subread-counts.txt $args->{bam}", $args->{logfile});

}

### Subreadplot ### - nobody calls this function
sub subreadplot {

    my ($samples, $args, $stats) = @_;
    print "Generating subread plots\n" if ($args->{verbose});

    open OFILE, ">subread-filelist.txt" or die "cannot open subread-filelist.tx\
t: $!\n";
    foreach $sample (sort keys %{$samples}) {
        print OFILE "$args->{scratchfolder}/rnaseqQC-ss/$sample/subread-counts.txt.summary\t$sample\n";
    }
    close OFILE;
    @result = `subreadplot.pl -f subread-filelist.txt`;

    # process the stdout, grabbing some stats
    $statzone = 0;
    foreach $result (@result) {
	$statzone = 0 if ($result =~ /^STATS End/);
	if ($statzone) {
	    chomp $result;
	    ($sample, $key, $stat) = split /\t/, $result;
	    $stats->{$sample}{subread}{$key} = $stat;
	}
	$statzone = 1 if ($result =~ /^STATS Start/);
    }
    `mv subread.png $args->{outputfolder}`;

    # add to report
    &h2("Subread Plots");
    &image("subread.png");
    &space;
}

######################## Cuffquant ##########################
sub cuffquant {
    print "Estimating abundances with cuffquant\n";

    my ($args, $stats) = @_;

    my $frags = ($stats->{subsample}{"#subsampledreads"} > 20000000) ? "--max-bundle-frags 10000000" : "";
    my $mask = ($args->{maskfile}) ? "--mask-file $args->{maskfile}" : "";

    my $cuffquantoptions = $args->{extraoptions}{cuffquant} // "--quiet --multi-read-correct $frags $mask";

    run("cuffquant $cuffquantoptions --no-update-check $args->{gtffile} $args->{bam} -o $args->{scratchfolder} -p $args->{threads}", $args->{logfile});

}

######################## GBS Trim ##########################
sub gbstrim {
    print "Running GBS Trim\n";

    my ($args, $stats) = @_;

    my $verbose = ($args->{verbose}) ? "--verbose" : "";
    my $padding = ($args->{gbspadding}) ? "--padding $args->{gbspadding}" : "";
    my $crop = ($args->{croplength}) ? "--minlength $args->{croplength} --maxlength $args->{croplength}" : "";
    $outputfastq = "$args->{scratchfolder}/R1.gbstrim.fastq";
    $logfile = "$args->{scratchfolder}/gbstrimR1.log";

    $resultr1 = &run("gbstrim.pl $verbose --enzyme $args->{enzyme} $padding --fastqfile $args->{R1} --outputfile $outputfastq $crop", $args->{logfile});

    open OFILE, ">$logfile" or die "cannot open gbstrim log file $logfile: $!\n";
    print OFILE $result;
    close OFILE;
    $args->{R1} = $outputfastq;

    if ($args->{R2}) {
	$outputfastq = "$args->{scratchfolder}/R2.gbstrim.fastq";
	$logfile = "$args->{scratchfolder}/gbstrimR2.log";
	$resultr2 = &run("gbstrim.pl $verbose --star --enzyme $args->{enzyme} $padding --fastqfile $args->{R2} --outputfile $outputfastq $crop", $args->{logfile});

	open OFILE, ">$logfile" or die "cannot open gbstrim log file $logfile: $!\n";
	print OFILE $resultr2;
	close OFILE;
	$args->{R2} = $outputfastq;
    }


    # process the stdout, grabbing some stats
    $statzone = 0;
    @results = split /\n/, $resultr1;
    foreach $result (@results) {
	$statzone = 0 if ($result =~ /^STATS End/);
	if ($statzone) {
	    chomp $result;
	    ($key, $stat) = split /\t/, $result;
	    $stats->{gbstrimR1}{$key} = $stat;
	}
	$statzone = 1 if ($result =~ /^STATS Start/);
    }

    # process the stdout, grabbing some stats
    if ($args->{R2}) {
	$statzone = 0;
	@results = split /\n/, $resultr2;
	foreach $result (@results) {
	    $statzone = 0 if ($result =~ /^STATS End/);
	    if ($statzone) {
		chomp $result;
		($key, $stat) = split /\t/, $result;
		$stats->{gbstrimR2}{$key} = $stat;
	    }
	    $statzone = 1 if ($result =~ /^STATS Start/);
	}
    }

    # resync files
    if ($args->{R2}) {
	$outputfastqr1 = "$args->{scratchfolder}/R1.gbstrim.resync.fastq";
	$outputfastqr2 = "$args->{scratchfolder}/R2.gbstrim.resync.fastq";
	&run("resync.pl $args->{R1} $args->{R2} $outputfastqr1 $outputfastqr2", $args->{logfile});
	$args->{R1} = $outputfastqr1;
	$args->{R2} = $outputfastqr2;
    }
	     
}

######################## all GBS Trim ##########################
sub allgbstrim {
    print "Running GBS Trim\n";
    
    push @{$args->{stageorder}}, "bgstrim";
    
    my ($samples, $args, $stats) = @_;

    my $myscratchfolder = "$args->{scratchfolder}/gbstrim";
    mkdir $myscratchfolder;

    my $verbose = ($args->{verbose}) ? "--verbose" : "";
    my $crop = ($args->{croplength}) ? "--minlength $args->{croplength} --maxlength $args->{croplength}" : "";
    open GP, "| parallel -j $args->{threads} >> $args->{logfile}" or die "Cannot open pipe to gnu parallel\n";
    foreach $sample (sort keys %{$samples}) {
	my $padding = "";
	$padding = "--padding $samples->{$sample}{padding}" if ($samples->{$sample}{padding} and ! $args->{headcrop});
	$outputfastq = "$myscratchfolder/${sample}_R1.trim.fastq";
#	print "enzyme: $args->{enzyme}\ninput: $samples->{$sample}{R1}{fastq}\noutput: $outputfastq\ncroplength $args->{croplength}\n";
	$command = "gbstrim.pl $verbose --enzyme $samples->{$sample}{enzyme} $padding --fastqfile $samples->{$sample}{R1}{fastq} --outputfile $outputfastq $crop &> $myscratchfolder/${sample}_R1.log";
#	print "$command\n";
	print GP "$command\n";
#       print $result;
	$samples->{$sample}{R1}{fastq} = $outputfastq;

	if ($samples->{$sample}{R2}{fastq}) {
	    my $padding2 = ($samples->{$sample}{padding2}) ? "--padding $samples->{$sample}{padding2}" : "";

	    $outputfastq = "$myscratchfolder/${sample}_R2.trim.fastq";
#	    print "enzyme: $args->{enzyme}\ninput: $samples->{$sample}{R1}{fastq}\noutput: $outputfastq\ncroplength $args->{croplength}\n";
	    $command = "gbstrim.pl $verbose --enzyme $samples->{$sample}{enzyme2} $padding2 --fastqfile $samples->{$sample}{R2}{fastq} --outputfile $outputfastq $crop &> $myscratchfolder/${sample}_R2.log";
#	    print "$command\n";
	    print GP "$command\n";
#           print $result;
	    $samples->{$sample}{R2}{fastq} = $outputfastq;
	}

    }
    close GP;

    # Read in and save stats
    foreach $sample (sort keys %{$samples}) {
	$ifile = "$myscratchfolder/${sample}_R1.log";
	if (!(open IFILE, $ifile)) {
	    print "Cannot open log file $ifile: $!\n";
	} else {
	    # process the stdout, grabbing some stats
	    $statzone = 0;
	    while ($result = <IFILE>) {
		$statzone = 0 if ($result =~ /^STATS End/);
		if ($statzone) {
		    chomp $result;
		    ($key, $stat) = split /\t/, $result;
		    $stats->{$sample}{gbstrimR1}{$key} = $stat;
		}
		$statzone = 1 if ($result =~ /^STATS Start/);
	    }
	}

	# handle R2 log, if paired-end data
	next unless ($samples->{$sample}{R2}{fastq});
	$ifile = "$myscratchfolder/${sample}_R2.log";
	if (!(open IFILE, $ifile)) {
	    print "Cannot open log file $ifile: $!\n";
	} else {
	    # process the stdout, grabbing some stats
	    $statzone = 0;
	    while ($result = <IFILE>) {
		$statzone = 0 if ($result =~ /^STATS End/);
		if ($statzone) {
		    chomp $result;
		    ($key, $stat) = split /\t/, $result;
		    $stats->{$sample}{gbstrimR2}{$key} = $stat;
		}
		$statzone = 1 if ($result =~ /^STATS Start/);
	    }
	}
    }

    # combine paired-end fastq files
    if ($args->{pe}) {
	foreach $sample (sort keys %{$samples}) {
	    $ofile = "$myscratchfolder/$sample.fastq";
	    `cat $samples->{$sample}{R1}{fastq} $samples->{$sample}{R2}{fastq} > $ofile`;
	    $samples->{$sample}{R1}{fastq} = $ofile;
	    $samples->{$sample}{R2}{fastq} = "";
	}
    } else {
	foreach $sample (sort keys %{$samples}) {
	    $ofile = "$myscratchfolder/$sample.fastq";
	    `mv $samples->{$sample}{R1}{fastq} $ofile`;
	    $samples->{$sample}{R1}{fastq} = $ofile;
	    $samples->{$sample}{R2}{fastq} = "";
	}	
    }
    $args->{fastqfolder} = "$myscratchfolder";
}

############################# Kraken #################################
# Requires: kraken
# Args: krakendb
sub kraken {
    print "Running Kraken\n";

    my ($args, $stats) = @_;
    my $myscratchfolder = "$args->{scratchfolder}/kraken";
    mkdir $myscratchfolder;

#    $result = &run("kraken --fastq-input --paired --preload --threads $args->{threads} --db $args->{krakendb} $args->{R1} $args->{R2} > $myscratchfolder/kraken.out", $args->{logfile});
    $result = &run("kraken --fastq-input --paired --preload --threads $args->{threads} --db $args->{krakendb} $args->{R1} $args->{R2} --output $myscratchfolder/kraken.out", $args->{logfile});

    $result = &run("kraken-report --db $args->{krakendb}  $myscratchfolder/kraken.out", $args->{logfile});
    # save stdout to a file
    $ofile = "$myscratchfolder/kraken.report";
    open OFILE, ">$ofile" or print "Cannot open kraken log $ofile: $!\n";
    print OFILE $result;
    close OFILE;

    # these reports can create very large files and I have not found a use for them yet.
    #kraken-translate --mpa-format --db $DBNAME ${SAMPLE_BASE[$i]}.kraken > ${SAMPLE_BASE[$i]}.labels
    #kraken-report --show-zeros --db $DBNAME ${SAMPLE_BASE[$i]}.kraken > ${SAMPLE_BASE[$i]}.reportzeros
    #kraken-mpa-report --db $DBNAME ${SAMPLE_BASE[$i]}.kraken > ${SAMPLE_BASE[$i]}

    # Generate some summaries of kraken output
    print "Generating Kraken summaries\n" if ($args->{verbose});
    my $originaldir = `pwd`;
    chdir $myscratchfolder;
    $result = &run("kraken_tabs.py $myscratchfolder", $args->{logfile});
    $result = &run("kraken_noindent.py $myscratchfolder", $args->{logfile});
    $result = &run("kraken_virus.py $myscratchfolder", $args->{logfile});
    $result = &run("kraken_virus_sort.py $myscratchfolder", $args->{logfile});
    chdir $originaldir;
}

### verify that a kraken database is at the specified location, make it's path absolute
sub krakendbcheck {
    my ($index) = @_;
    my $indexfile = "$index/database.idx";
    if (-e $indexfile) {
	$indexfile = abs_path($indexfile);
    } else {
	die "Unable to find KrakenDB $index";
    }
    ($name, $path) = fileparse($indexfile);
    return "$path/";
}


############################## Kraken Plot #######################
sub krakenplot {
    my ($samples, $args, $stats) = @_;

    print "\nGenerating Kraken Plots\n" if ($args->{verbose});

# Create Kraken report_summary file
# ---------------------------------------------------------------------
# Extract relevant percent and count data from the kraken report file

    $myscratchfolder = "$args->{scratchfolder}/kraken";
    mkdir $myscratchfolder;
    $krakenreport = "$myscratchfolder/kraken.report_summary";
    foreach $sample (sort keys %{$samples}) {
	`echo "$sample" >> $krakenreport`;

	# Unclassified Reads (taxonID: 0)
	$stats->{$sample}{kraken}{percentunclassified} = 0;
	$stats->{$sample}{kraken}{numberunclassified} = 0;
	`awk '\$5 == "0"' $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.report >> $krakenreport`;
	$result = `awk '\$5 == "0"' $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.report`;
	if ($result) {
	    ($percent, $number, $numberdirect, $rank, $taxonid) = split ' ', $result;
	}

	# Unassigned Reads (taxonID: 1)
	$stats->{$sample}{kraken}{totalinputreads} = 0;
	$stats->{$sample}{kraken}{percentunassignedroot} = 0;
	$stats->{$sample}{kraken}{numberunassignedroot} = 0;
	$result = `awk '\$5 == "1"' $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.report`;
	if ($result) {
	    ($percent, $number, $numberdirect, $rank, $taxonid) = split ' ', $result;
	    $stats->{$sample}{kraken}{totalinputreads} = $number + $stats->{$sample}{kraken}{numberunclassified};
	    $stats->{$sample}{kraken}{percentunassignedroot} = round100($numberdirect / $stats->{$sample}{kraken}{totalinputreads} * 100);
	    $stats->{$sample}{kraken}{numberunassignedroot} = $numberdirect;
	}

	# Cellular organisms (taxonID: 1131567)
	$stats->{$sample}{kraken}{percentcellular} = 0;
	$stats->{$sample}{kraken}{numbercellular} = 0;
	$result = `awk '\$5 == "131567"' $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.report`;
	if ($result) {
	    ($percent, $number, $numberdirect, $rank, $taxonid) = split ' ', $result;
	    $stats->{$sample}{kraken}{percentcellular} = round100($numberdirect / $stats->{$sample}{kraken}{totalinputreads} * 100);
	    $stats->{$sample}{kraken}{numbercellular} = $numberdirect;
	}

	# Eukaryota (taxonID: 2759)
	$stats->{$sample}{kraken}{percenteukaryta} = 0;
	$stats->{$sample}{kraken}{numbereukaryota} = 0;
	$result = `awk '\$5 == "2759"' $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.report`;
	if ($result) {
	    ($percent, $number, $numberdirect, $rank, $taxonid) = split ' ', $result;
	    $stats->{$sample}{kraken}{percenteukaryta} = round100($numberdirect / $stats->{$sample}{kraken}{totalinputreads} * 100);
	    $stats->{$sample}{kraken}{numbereukaryota} = $numberdirect;
	}

	# Host Reads (pig == Sus scrofa == taxonID: 9823)
	$stats->{$sample}{kraken}{percenthost} = 0;
	$stats->{$sample}{kraken}{numberhost} = 0;
	`awk '\$5 == "$args->{taxonid}"' $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.report >> $krakenreport`;
	$result = `awk '\$5 == "$args->{taxonid}"' $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.report`;
	if ($result) {
	    ($percent, $number, $numberdirect, $rank, $taxonid) = split ' ', $result;
	    $stats->{$sample}{kraken}{percenthost} = $percent;
	    $stats->{$sample}{kraken}{numberhost} = $number;
	}
	
	# Bacteria Reads (taxonID: 2)
	$stats->{$sample}{kraken}{percentbacteria} = 0;
	$stats->{$sample}{kraken}{numberbacteria} = 0;
	`awk '\$5 == "2"' $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.report >> $krakenreport`;
	$result = `awk '\$5 == "2"' $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.report`;
	if ($result) {
	    ($percent, $number, $numberdirect, $rank, $taxonid) = split ' ', $result;
	    $stats->{$sample}{kraken}{percentbacteria} = $percent;
	    $stats->{$sample}{kraken}{numberbacteria} = $number;
	}

	# Viruses Reads (taxonID: 10239)
	$stats->{$sample}{kraken}{percentvirus} = 0;
	$stats->{$sample}{kraken}{numbervirus} = 0;
	`awk '\$5 == "10239"' $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.report >> $krakenreport`;
	$result = `awk '\$5 == "10239"' $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.report`;
	if ($result) {
	    ($percent, $number, $numberdirect, $rank, $taxonid) = split ' ', $result;
	    $stats->{$sample}{kraken}{percentvirus} = $percent;
	    $stats->{$sample}{kraken}{numbervirus} = $number;
	}
    }

# ---------------------------------------------------------------------
# Create a file to be used by Krona (pie charts)
# ---------------------------------------------------------------------

# Extract taxonomy number from kraken results
# Parse the Kraken output file (keeping only the queryID, NCBI taxonomyID columns)
    # Columns 2 and 3 of the original Kraken results are kept and saved to new file
# Also, create a space-delimited string of all the "kraken_for_krona" filenames
    $files_for_krona = "";
    foreach $sample (sort keys %{$samples}) {
	`cut -f2,3 $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.out > $myscratchfolder/$sample.kraken_for_krona`;
#	`tail -n +2 $args->{scratchfolder}/singlesamples/$sample/kraken/kraken.out | cut -f2,3 > $myscratchfolder/$sample.kraken_for_krona`;
	$files_for_krona .= "$myscratchfolder/$sample.kraken_for_krona ";
    }

# ---------------------------------------------------------------------
# Run Krona (make pie charts)
# ---------------------------------------------------------------------

    # make a krona dir 
    mkdir "$args->{scratchfolder}/krona";

# Run the Krona script for importing taxonomy data (for all files)
# Place output files into a different directory (../krona) because pwd == kraken.
    print "files for krona: $files_for_krona\n";
    $result = run("ktImportTaxonomy -i -k -o $args->{scratchfolder}/krona/krona.html $files_for_krona", $args->{logfile});
    print $result;

    # move and rename plots
    `cp $args->{scratchfolder}/krona/krona.html $args->{outputfolder}`;
    
    # add to report
    &h1("Kraken");
    &h1("Krona");
#    &image("${name}_heatmapplot.png");
    &line("`Krona Plot <Resources/krona.html>`_"); 

    &h2("Kraken Reports"); 
    foreach $sample (sort keys %{$samples}) {
	print $reportFH "- $sample `No-indent <Resources/kraken/$sample" . ".report.noindent.txt>`_";
	print $reportFH " `Tabs <Resources/kraken/$sample" . ".report.tabs.txt>`_";
	print $reportFH " `Virus <Resources/kraken/$sample" . ".report.virus.sort.txt>`_";
	print $reportFH "\n";
    }
    &space;

}

############################# Gunzip #################################
sub allgunzip {
    my ($samples, $args, $stats) = @_;

    foreach $sample (keys %{$samples}) {
	$args->{gz} = 1 if ($samples->{$sample}{R1}{fastq} =~ /\.gz$/);
	$args->{gz} = 1 if (($samples->{$sample}{R2}{fastq}) and ($samples->{$sample}{R2}{fastq} =~ /\.gz$/));
    }
    if (! $args->{gz}) {
	print "Skipping gunzip, files not compressed\n" if ($args->{verbose});
	return;
    }

    print "Running gunzip on all samples\n";

    push @{$args->{stageorder}}, "gunzip";

    # gunzip files
    my $myscratchfolder = "$args->{scratchfolder}/fastq-gunzip";
    mkdir $myscratchfolder;
    $gpout = "$myscratchfolder/gp.out";
    $progress = $args->{verbose} ? "--progress" : "";
    open GP, "| parallel $progress -j $args->{threads} > $gpout" or die "Cannot open pipe to gnu parallel\n";
#    open GP, ">$myscratchfolder/parallel.test" or die "Cannot open pipe to gnu parallel\n";

    # send commands to gnu parallel
    foreach $sample (keys %{$samples}) {
	my ($fname, $path) = fileparse($samples->{$sample}{R1}{fastq}, (".gz"));
	$newfile = "$myscratchfolder/$fname";
	print GP "gunzip -c $samples->{$sample}{R1}{fastq} > $newfile\n";
	$samples->{$sample}{R1}{fastq} = $newfile;
	if (defined($samples->{$sample}{R2}{fastq})) {
	    ($fname, $path) = fileparse($samples->{$sample}{R2}{fastq}, (".gz"));
	    $newfile = "$myscratchfolder/$fname";
	    print GP "gunzip -c $samples->{$sample}{R2}{fastq} > $newfile\n";
	    $samples->{$sample}{R2}{fastq} = $newfile;
	}
    }

    close GP;
    print "error code: $?\n" if ($?);
    $args->{fastqfolder} = "$myscratchfolder";

}


############### Single sample setup and helper functions ##############

### Set up a scratch folder, add symlinks to fastq files
sub scratchfoldersetup {
    my ($args) = @_;

    my $scratchfolder = $args->{scratchfolder};

    ### set up scratch folders
    if (-e $scratchfolder) {
	print "scratch folder $scratchfolder exists, deleting it\n";
        `rm -r $scratchfolder`;
	die "Error: Cannot remove existing scratchfolder $scratchfolder\n" if (-e $scratchfolder);
    }
    mkdir $scratchfolder or die "Cannot create scratchfolder $scratchfolder: $!\n";
    chdir $scratchfolder;

    # make local copy of fastq files
    if ($args->{R1} =~ /\.gz$/) {
	# uncompress gzipped fastq files
	run("gunzip -c $args->{R1} > $scratchfolder/R1.fastq", $args->{logfile});
	$args->{R1} = "$scratchfolder/R1.fastq";
    } else {
	# create fastq file symlinks
	`ln -s $args->{R1} $scratchfolder/R1.fastq`;
	$args->{R1} = "$scratchfolder/R1.fastq";
    }
    if ($args->{pe}) {

	if ($args->{R2} =~ /\.gz$/) {
	    # uncompress gzipped fastq files
	    run("gunzip -c $args->{R2} > $scratchfolder/R2.fastq", $args->{logfile});
	    $args->{R2} = "$scratchfolder/R2.fastq";
	} else {
	    # create fastq file symlinks
	    `ln -s $args->{R2} $scratchfolder/R2.fastq`;
	    $args->{R2} = "$scratchfolder/R2.fastq";
	}
    }
}

# write the %stats hash out to an easily-parsed text file
sub writestats {
    print "Writing stats\n";

    my ($args, $stats) = @_;

    my $ofile = "$args->{scratchfolder}/stats.txt";
    open OFILE, ">$ofile" or die "Cannot open stats file $ofile\n";
    foreach $stage (keys %{$stats}) {
	foreach $stat (keys %{$stats->{$stage}}) {
	    print OFILE "$stage\t$stat\t$stats->{$stage}{$stat}\n";
	}
    }
    close OFILE;
}

############################################################################
############################# All Sample Scripts ###########################
############################################################################

### bam2fastq ###
sub bam2fastq {
    my ($samples, $args, $stats, $vcf) = @_;

    print "Running bam2fastq\n" if ($args->{verbose});
    $myscratchfolder = "$args->{scratchfolder}/bam2fastq";
    mkdir $myscratchfolder;
    mkdir "$myscratchfolder/aligned";
    mkdir "$myscratchfolder/unaligned";    

    $gpout = "$myscratchfolder/bam2fastq.gpout";
    open GP, "| parallel -j $args->{threadspersample} > $gpout" or die "Cannot open pipe to gnu parallel\n"; 
    open FILE, ">$myscratchfolder/commands.txt" or die "Cannot open pipe to gnu parallel\n";

    # TODO: add -f or -F flags to include only primary alignments, or only unaligned reads
    $filter = "-f 4";
    # then add this to align-pipeline
    foreach $sample (keys %{$samples}) {
	# unaligned reads
	$command = "samtools fastq -f 4 -1 $myscratchfolder/unaligned/${sample}_R1.fastq.gz -2 $myscratchfolder/unaligned/${sample}_R2.fastq.gz $args->{scratchfolder}/bam/$sample.bam \n";
	print GP $command;
	print FILE $command;
	# aligned reads: don't print if unmapped or not primary alignment
	$command = "samtools fastq -F 4 -F 256 -1 $myscratchfolder/aligned/${sample}_R1.fastq.gz -2 $myscratchfolder/aligned/${sample}_R2.fastq.gz $args->{scratchfolder}/bam/$sample.bam \n";
	print GP $command;
	print FILE $command;
    }
    close GP;
    close IFILE;

}


# run freebayes on all samples (bam files)
sub freebayes {
    print "Calling variants with Freebayes\n";

    my ($samples, $args, $stats) = @_;

    push @{$args->{stageorder}}, "freebayes";

    my $myscratchfolder = "$args->{scratchfolder}/freebayes";
    mkdir $myscratchfolder;
    my $vcffolder = "$myscratchfolder/vcffiles";
    mkdir $vcffolder;
    $vcfrawfile = "$myscratchfolder/variants.raw.vcf";
    $vcffile = "$myscratchfolder/variants.vcf";
    $vcffilterfile = "$myscratchfolder/variants.filt.vcf";
    my $ploidy = $args->{ploidy} ? "--ploidy $args->{ploidy}" : "";

    $bamlist = "";
    # make local bam folder to keep freebayes command line length within unix limits
    mkdir "$myscratchfolder/bams";
    foreach $sample (keys %{$samples}) {
#	print "Freebayes bam $sample $stats->{$sample}{files}{bam}\n";
	`ln -s $stats->{$sample}{files}{bam} $myscratchfolder/bams/$sample.bam`;
	`ln -s $stats->{$sample}{files}{bam}.bai $myscratchfolder/bams/$sample.bam.bai`;
	$bamlist .= " bams/$sample.bam";
#	$bamlist .= " $stats->{$sample}{files}{bam}";
    }
#    print "BAMlist: $bamlist end\n";

    my $originaldir = `pwd`;
    chdir $myscratchfolder;

    $freebayesthreads = $args->{freebayesthreads} // $args->{threads};

    if ($freebayesthreads > 1) {

	$gpout = "$myscratchfolder/freebayes.gpout";
	open GP, "| parallel -j $freebayesthreads > $gpout" or die "Cannot open pipe to gnu parallel\n"; 
	open FILE, ">$myscratchfolder/commands.txt" or die "Cannot open pipe to gnu parallel\n";
	$myfai = "$args->{referencefasta}.fai";	
	open IFILE, $myfai or die "Cannot open genome fasta index $myfai: $!\n";
	$scaffoldcount = 0;
	while ($line = <IFILE>) {
	    chomp $line;
	    my ($scaffold, $length) = split ' ', $line;
#	    $region = "$scaffold:0-$length";
	    $region = "$scaffold";
	    $nicescaffold = $scaffold;
	    $nicescaffold =~ s/\W/./g;
	    $scaffoldcount++;
	    if ($args->{resume} && -e "$vcffolder/$nicescaffold.complete") {
		print "Skipping completed scaffold $nicescaffold\n";
		next;
	    }
	    $output = "$vcffolder/$nicescaffold.vcf";
	    $command = "freebayes --genotype-qualities $ploidy --fasta-reference $args->{referencefasta} $bamlist --region \"$region\" > $output && touch $vcffolder/$nicescaffold.complete\n";
	    print GP $command;
	    print FILE $command;
	}
	close GP;
	close IFILE;

	# check that all freebayes completed successfully
	$completecount = `ls $vcffolder/*.complete | wc -l`;
	chomp $completecount;
	die "Freebayes error, freebayes was successful on only $completecount of $scaffoldcount scaffolds\n" if ($completecount != $scaffoldcount);

	# combine all vcf files into one
	&run("vcf-merge.pl $vcffolder $vcfrawfile", $args->{logfile});
#	$result = `gzip $tmpvcf`;
	print $result if ($args->{verbose});
	
    } else {
	&run("freebayes --genotype-qualities $ploidy --fasta-reference $args->{referencefasta} $bamlist > $vcfrawfile", $args->{logfile});
    }
    chdir $originaldir;

    # basic filter: variants with QUAL > 20
    &run("vcffilter -f \"QUAL > 20\" $vcfrawfile > $vcffile", $args->{logfile});

    # stringent filter: remove samples with > 5% missing genotypes, and variants with missing calls in more than 5% of samples
    $samplecutoff = $args->{samplecutoff} ? "--samplecutoff $args->{samplecutoff}" : "";
    $markercutoff  = $args->{markercutoff} ? "--markercutoff $args->{markercutoff}" : "";
    $maf  = $args->{maf} ? "--maf $args->{maf}" : "";
    &run("vcffilter.pl $samplecutoff $markercutoff $maf --vcffile $vcffile --output $vcffilterfile", $args->{logfile});

    $args->{vcffilter} = $vcffilterfile;
    $args->{vcf} = $vcffile;

}

### freebayes Plot ###
sub freebayesplot {
    my ($samples, $args, $stats) = @_;

# the freebayes subroutine already generates the plots, here we just add the 
# plots to the report, this provides control over where in the report
# the plots appear

    # add to report
    &h2("Freebayes Plots");  
    &image("pcaplot.png");
    &image("pcaplot-filt.png");
    &space;

}


### vcfMetrics Plot ###
sub vcfmetricsplot {
    my ($samples, $args, $stats, $vcf) = @_;

    print "Running vcfmetricsplot\n" if ($args->{verbose});
    my $vcffile = $vcf // $args->{vcf};

    ($prefix) = fileparse($vcffile);
    ($prefix) = split /vcf$/, $prefix;

    @result = run("vcfmetricsplot.pl --vcffile $vcffile --prefix $prefix", $args->{logfile});

    # process the stdout, grabbing some stats
    $statzone = 0;
    print "Looking for stats...\n";
    print @result;
    foreach $result (@result) {
        $statzone = 0 if ($result =~ /^STATS End/);
        if ($statzone) {
            chomp $result;
	    print "result: $result\n";
            ($sample, $key, $stat) = split /\t/, $result;
	    $key = $prefix . "." . $key;
            $stats->{$sample}{vcfmetrics}{$key} = $stat;
        }
        $statzone = 1 if ($result =~ /^STATS Start/);
    }

    my $firstsample = (keys %{$stats})[0];
    ($vcfname) = fileparse($vcffile);

    # generate chrom plot
    if ($args->{referencefasta}) {
	run("chromplot.pl --vcffile $vcffile --referencefai $args->{referencefasta}.fai --prefix $prefix", $args->{logfile});
    }

    # copy all plots to outputfolder
    `mv *.png $args->{outputfolder}`;
    `mv *.space $args->{outputfolder}`;

    # create data file of marker numbers at different filter cutoffs
    $ofile = $prefix . "markers";
    if ($args->{referencefasta}) {
	$thin = 100;
    } else {
	$thin = 0
    }
    open OFILE, ">$ofile";
    foreach $cutoff (.99, .95, .9, .8, .7, .6, .5) {
	$result = `vcffilter.pl --vcffile $vcffile --outputfile tmp.vcf --samplecutoff $cutoff --markercutoff $cutoff --thin $thin 2>/dev/null`;
	$result =~ /kept (\d+) out of/;
	$markers = $1;
	if ($thin == 0) {
	    $markers = `grep -v "^#" tmp.vcf | cut -f 3 | cut -f 1 -d '_' | sort | uniq | wc -l`;
	    chomp $markers;
	}
	$result =~ /(\d+) bad samples filtered out from a total of (\d+) samples/;
	$bad = $1;
	$total = $2;
	$good = $total - $bad;
	print OFILE "$cutoff\t$markers\t$good\n";
    }
    close OFILE;
    `mv *.markers $args->{outputfolder}`;

    # add to report
    &h1("VCF Metrics Plots: $vcfname");
#    &fieldlist("Bait territory", "$stats->{$firstsample}{hsmetrics}{BAIT_TERRITORY}bp");
#    &fieldlist("Target territory", "$stats->{$firstsample}{hsmetrics}{TARGET_TERRITORY}bp");
#    &fieldlist("Genome size", "$stats->{$firstsample}{hsmetrics}{GENOME_SIZE}bp");
    &space;
    &image($prefix . "pca.png");
    &image($prefix . "sample-depth.png");
    &image($prefix . "sample-miss.png");
#    &image($prefix . "sample-het.png");
    &image($prefix . "site-depth.png");
    &image($prefix . "site-miss.png");
    &image($prefix . "site-space.png") if ($args->{referencefasta} or $args->{bwaindex});
    &image($prefix . "chromplot.png") if ($args->{referencefasta});
    &image($prefix . "chromcountplot.png") if ($args->{referencefasta});
#    &image($prefix . "site-qual.png"); # this isn't very helpful
    &space;
}

# set up scratch and output folder
# read in samplesheet, make sure it is properly formed
sub samplesheetandfoldersetup {
    print "Processing samplesheet and setting up folders\n";

    my ($samples, $args, $stats, $pipeline) = @_;

    push @{$args->{stageorder}}, "general";

    ### set up folders
    if ($args->{samplesheet}) {
	die "Cannot find samplesheet $args->{samplesheet}\n" if (! -e $args->{samplesheet});
	$args->{samplesheet} = abs_path($args->{samplesheet});
    }
#    if ($pipeline ne "umgc") { # set up fastq folder
	die "Fastq folder $args->{fastqfolder} not found\n" if (!-e $args->{fastqfolder});
	$args->{fastqfolder} = abs_path($args->{fastqfolder});
#    }
	my ($runname) = ($args->{fastqfolder} =~ /.*\/([^\/]+)\/?/);
	$args->{runname} = $runname unless ($args->{runname});

    my $user = `whoami`;
    chomp $user;

    $args->{scratchfolder} = $args->{scratchfolder} // "$SCRATCHGLOBAL/$user-pipelines/$pipeline-$args->{runname}";

    if (!$args->{resume}) {
	if (-e $args->{scratchfolder}) {
	    print "scratch folder $args->{scratchfolder} exists\n";
	    $|++;
	    for $i (4..5) { 
		$remaining = 5-$i;
		print "\rDELETING $args->{scratchfolder} in $remaining seconds";
		sleep 1;
	    }
	    print "\n";
	    `rm -r $args->{scratchfolder}`;
	    
	    die "Error: Cannot remove existing scratchfolder $args->{scratchfolder}\n" if (-e $args->{scratchfolder});
	}
	
	`mkdir --parents $args->{scratchfolder}`;
	die "Cannot create temporary folder $args->{scratchfolder}\n" if (! -e $args->{scratchfolder});
	$args->{scratchfolder} = abs_path($args->{scratchfolder});
	mkdir "$args->{scratchfolder}/allsamples";
	
	$args->{outputfolder} = "$args->{scratchfolder}/output";
	mkdir "$args->{outputfolder}";
	die "Cannot create output folder $args->{outputfolder}\n" if (! -e $args->{outputfolder});

	mkdir "$args->{scratchfolder}/fastq";
	$args->{logfile} = "$args->{scratchfolder}/log";
	
	### Process samplesheet
	# create samplesheet if not provided
	if (!defined($args->{samplesheet})) {
	    print "No samplesheet provided, attempting to create one\n";
	    $args->{samplesheet} = "$args->{scratchfolder}/samplesheet.txt";
	    if ( ($pipeline eq "umgc") or ($pipeline eq "illumina") ) { # require illumina fastq filenames
		$result = `createsamplesheet.pl -i -f $args->{fastqfolder} -o $args->{samplesheet}`;
	    } elsif ($pipeline eq "metagenomics") {
		if ($args->{emp} or $args->{variableregion}) {
		    $result = &run("createsamplesheet.pl -p $ENV{GOPHER_PIPELINES}/resources/16s-primers.fa -q -e -f $args->{fastqfolder} -o $args->{samplesheet}", $args->{logfile});
		} else {
		    $result = &run("createsamplesheet.pl -p $ENV{GOPHER_PIPELINES}/resources/16s-primers.fa -q -f $args->{fastqfolder} -o $args->{samplesheet}", $args->{logfile});
		}
	    } else {
		$result = `createsamplesheet.pl -f $args->{fastqfolder} -o $args->{samplesheet}`;
	    }
	    print $result;
	    die "Samplesheet could not be created\n" if (! -e $args->{samplesheet});
	} else { # copy over the specified samplesheet
	    $result = `cp $args->{samplesheet} $args->{scratchfolder}/samplesheet.txt`;
	    print $result;
	    $args->{samplesheet} = "$args->{scratchfolder}/samplesheet.txt";
	}
    } else { # resumeing a run
	die "Attempting to resume run but scratch folder doesn't exist: $args->{scratchfolder}\n(Re-run without the --resume option)\n" unless (-e $args->{scratchfolder});
	# copy over samplesheet if it is newer than the copy in scratch
	if (defined($args->{samplesheet})) {
	    if (-M $args->{samplesheet} < -M "$args->{scratchfolder}/samplesheet.txt") {
		$result = `cp $args->{samplesheet} $args->{scratchfolder}/samplesheet.txt`;
		print $result;
	    }
	}
	$args->{outputfolder} = $args->{outputfolder} // "$args->{scratchfolder}/output";
	$args->{outputfolder} = abs_path($args->{outputfolder}) if ($args->{outputfolder});
	$args->{samplesheet} = "$args->{scratchfolder}/samplesheet.txt";
	$args->{outputfolder} = "$args->{outputfolder}";
    }	

    chdir $args->{scratchfolder};

    fixsamplesheet($args->{samplesheet});

    # read in the samplesheet
    open MFILE, $args->{samplesheet} or die "Cannot open samplesheet " . $args->{samplesheet};
    $header = <MFILE>;
    chomp $header;
    @header = split /\t/, lc($header);
    if ($header =~ /fastqR2/i) {
	$args->{pe} = 1;
	print "Paired-end read dataset\n";
    } else {
	$args->{pe} = 0;
	print "Single-end read dataset\n";
    }
    while ($line = <MFILE>) {
	chomp $line;
	@line = split /\t/, $line;
	$sample = $line[0];
	for $i (1..$#line) {
	    $samplesheetdata{$sample}{$header[$i]} = $line[$i] // "";
	}
    }
    close MFILE;

    # verify samplesheet file has required columns
    die "Samplesheet must contain a \"fastqR1\" column\n" if ($header !~ /fastqR1/i);
    die "Samplesheet must contain a \"fastqR2\" column\n" if (($header !~ /fastqR2/i) && ($args->{pe}));
    if ($pipeline eq "metagenomics") {
	die "Samplesheet (mapping file) must contain a \"ReversePrimer\" column\n" if (($header !~ /ReversePrimer/) && ($args->{pe}));
    }

    # print out samplesheet info for debugging
    if (0) { 
    foreach $sample (keys %samplesheetdata) {
	print "$sample";
	foreach $key (sort keys %{$samplesheetdata{$sample}}) {
	    print ":$key=$samplesheetdata{$sample}{$key}";
	}
	print "\n";
    }
    }

    ### For each sample save info from the samplesheet for later use
    foreach $sample (keys %samplesheetdata) {
	$file = $samplesheetdata{$sample}{fastqr1};
	if (! -e "$args->{fastqfolder}/$file") {
	    print "Could not find file $args->{fastqfolder}/$file\n";
	    next;
	}
	if (-z "$args->{fastqfolder}/$file") {
	    print "Excluding sample $sample from analysis due to empty fastq file $file\n";
	    next;
	}
	if ($file =~ /\.gz$/) {
	    $result = `gunzip -c $args->{fastqfolder}/$file | head -n 4`;
	    $result =~ s/\s+//g;
	    if ($result eq "") {
		print "Excluding sample $sample from analysis due to empty fastq file $file\n";
		next;
	    }
	}
	$fullname = $sample . "_R1.fastq";
	$fullname .= ".gz" if ($file =~ /\.gz$/);
	$samples->{$sample}{R1}{fastq} = "$args->{scratchfolder}/fastq/$fullname";
	push @symlinks, " ln -s $args->{fastqfolder}/$file $args->{scratchfolder}/fastq/$fullname" unless (-e "$args->{scratchfolder}/fastq/$fullname" and $args->{resume});
	if ($args->{pe}) {
	    $file = $samplesheetdata{$sample}{fastqr2};
	    if (! -e "$args->{fastqfolder}/$file") {
		print "Could not find file $args->{fastqfolder}/$file\n";
	    }
	    if (-z "$args->{fastqfolder}/$file") {
		print "Excluding sample $sample from analysis due to empty fastq file $file\n";
		next;
	    }
	    if ($file =~ /\.gz$/) {
		$result = `gunzip -c $args->{fastqfolder}/$file | head -n 4`;
		$result =~ s/\s+//g;
		if ($result eq "") {
		    print "Excluding sample $sample from analysis due to empty fastq file $file\n";
		    next;
		}
	    }
	    $fullname = $sample . "_R2.fastq";
	    $fullname .= ".gz" if ($file =~ /\.gz$/);
	    $samples->{$sample}{R2}{fastq} = "$args->{scratchfolder}/fastq/$fullname";
	    push @symlinks, " ln -s $args->{fastqfolder}/$file $args->{scratchfolder}/fastq/$fullname" unless (-e "$args->{scratchfolder}/fastq/$fullname" and $args->{resume});
	}

	if ($pipeline eq "metagenomics") {
	    $samples->{$sample}{forwardprimer} = $samplesheetdata{$sample}{linkerprimersequence} // "";
	    $samples->{$sample}{reverseprimer} = $samplesheetdata{$sample}{reverseprimer} // "";
	} 

	if (($pipeline eq "gbs") or ($pipeline eq "gbsref") or ($pipeline eq "tassel") or ($pipeline eq "stacks2")) {
	    $samples->{$sample}{enzyme} = $args->{enzyme} // $samplesheetdata{$sample}{enzyme} // "";
	    $samples->{$sample}{enzyme} = lc($samples->{$sample}{enzyme});
	    $samples->{$sample}{enzyme2} = $args->{enzyme2} // $samplesheetdata{$sample}{enzyme2} || $samples->{$sample}{enzyme} // "";
	    $samples->{$sample}{enzyme2} = lc($samples->{$sample}{enzyme2});
	    unless ($args->{nopadding}) {
		$samples->{$sample}{padding} = $args->{padding} // $samplesheetdata{$sample}{padding} || $padding{$samples->{$sample}{enzyme} . "-r1"} // "";
		$samples->{$sample}{padding2} = $args->{padding2} // $samplesheetdata{$sample}{padding2} || $padding{$samples->{$sample}{enzyme2} . "-r2"} // $samples->{$sample}{padding} // "";
	    }
	} 

   }

    foreach $symlink (@symlinks) {
	`$symlink`;
    }

    $args->{samplenumber} = scalar (keys %{$samples});
    print "Samples found: $args->{samplenumber}\n";
    die "No fastq files could be found\n" if ($args->{samplenumber} == 0);
    foreach $sample (keys %{$samples}) {
	$stats->{$sample}{general}{samplename} = $sample;
	$stats->{$sample}{general}{runname} = $args->{runname};
	$readlength = `fastqreadlength.pl $samples->{$sample}{R1}{fastq}`;
	chomp $readlength;
	$stats->{$sample}{general}{readlength} = $readlength;
    }
    
}

######################## Start Report ##########################
sub startreport {
    print "\nStarting Analysis Report\n";

    my ($samples, $args, $stats, $title, $brand) = @_;

    $title = "Analysis Report" unless (defined($title));
    $brand = "UMII" unless (defined($brand));

    open $reportFH, ">$args->{outputfolder}/AnalysisReport.rst" or die "Cannot open file $args->{outputfolder}/AnalysisReport.rst for writing: $!\n";

    if ($brand eq "UMGC") {
	# add UMGC logo?
    } else {
	&image("UMII-full-banner.png");
	&space;
    }
    &title($title);
#    &subtitle("Report");
    &fieldlist("Analysis run name", $args->{runname});
    &fieldlist("Date generated", $args->{startdate});
    &fieldlist("Samples", $args->{samplenumber});
    if ($args->{pe}) {
	&fieldlist("Library type", "Paired-end");
    } else {
	&fieldlist("Library type", "Single-end");
    }
    &fieldlist("Fastq folder", $args->{fastqfolder});
    &fieldlist("Commandline", "$0 $args->{options}");
    &fieldlist("Program versions", "`programversions.txt <Resources/programversions.txt>`_");
    &fieldlist("Sample metrics", "`metrics.html <Resources/metrics.html>`_ `metrics.txt <Resources/metrics.txt>`_");
    &space;

#    close $reportFH;
}

############################# Get singlesample stats ###################
sub getsinglesamplestats {
    print "\nGetting singlesample stats\n";
    my ($samples, $args, $stats) = @_;

    # Slurp up fastqQC stats ###
    foreach $sample (keys %{$samples}) {
	$ifile = "$args->{scratchfolder}/singlesamples/$sample/stats.txt";
	open IFILE, $ifile or warn "Cannot open stats file $ifile: $!\n";
	while ($line = <IFILE>) {
	    chomp $line;
	    ($stage, $stat, $value) = split /\t/, $line;
	    $stats->{$sample}{$stage}{$stat} = $value;
	}
    }

}

################### Compile Folder #######################
sub compilefolder {
    my ($samples, $args, $stats, $folder, $sourcefile, $destext) = @_;
    
    print "Compiling $folder folder\n" if ($args->{verbose});

    my $myscratchfolder = "$args->{scratchfolder}/$folder";
    mkdir $myscratchfolder;
    foreach $sample (sort keys %{$samples}) {
#	my $source = "$args->{scratchfolder}/singlesamples/$sample/$sourcefile"; # absolute paths
	my $source = "../singlesamples/$sample/$sourcefile"; # relative paths
	if (! -e $source) {
	    print "Cannot find file $source\n";
	    next;
	}
#	my $dest = "$myscratchfolder/$sample$destext"; # absolute paths
	my $dest = "../$folder/$sample$destext"; # relative paths
#	`mv --no-clobber $source $dest 2>&1`; # actually move the files
	$result = `ln -s $source $dest 2>&1`; # just create symlinks to them
	print $result;
    }
}

################### fastqcplot #######################
sub fastqcplot {
    print "\nGenerating FastQC Plot\n";

    my ($samples, $args, $stats) = @_;

    # move FastQC zips and htmls to fastqc folder
    my $ss = ($args->{subsample}) ? ".subsample" : ""; 
    my $ps = ($args->{prinseq}) ? ".prinseq" : ""; 
    $source = "fastqc/R1" . $ss . "_fastqc.html";
    &compilefolder($samples, $args, $stats, "fastqc", $source, "_R1.html");
    $source = "fastqc/R1" . $ss . "_fastqc.zip";
    &compilefolder($samples, $args, $stats, "fastqc", $source, "_R1.zip");
    if ($args->{pe}) {
	$source = "fastqc/R2" . $ss . "_fastqc.html";
	&compilefolder($samples, $args, $stats, "fastqc", $source, "_R2.html");
	$source = "fastqc/R2" . $ss . "_fastqc.zip";
	&compilefolder($samples, $args, $stats, "fastqc", $source, "_R2.zip");
    }
    # copy the html files to the output folder
    `cp -rL $args->{scratchfolder}/fastqc $args->{outputfolder}`;
    `rm $args->{outputfolder}/fastqc/*.zip`;

    ### Create filtered FastQC folder ###
    if ($args->{qualitycontrol}) {
	my $ss = ($args->{subsample}) ? ".subsample" : ""; 
	$source = "fastqc-trim/R1" . $ss . $ps . ".trim_fastqc.html";
	&compilefolder($samples, $args, $stats, "fastqc_filtered", $source, "_R1.html");
	$source = "fastqc-trim/R1" . $ss . $ps . ".trim_fastqc.zip";
	&compilefolder($samples, $args, $stats, "fastqc_filtered", $source, "_R1.zip");
	if ($args->{pe}) {
	    $source = "fastqc-trim/R2" . $ss . $ps . ".trim_fastqc.html";
	    &compilefolder($samples, $args, $stats, "fastqc_filtered", $source, "_R2.html");
	    $source = "fastqc-trim/R2" . $ss . $ps . ".trim_fastqc.zip";
	    &compilefolder($samples, $args, $stats, "fastqc_filtered", $source, "_R2.zip");
	}
	# copy the html files to the output folder
	`cp -rL $args->{scratchfolder}/fastqc_filtered $args->{outputfolder}`;
	`rm $args->{outputfolder}/fastqc_filtered/*.zip`;
    }

    ### fastqqualityplot ###
    print "Running fastqqualityplot.pl\n" if ($args->{verbose});
    open OFILE, ">fastqqualityplot-filelist.txt" or die "cannot open fastqqualityplot-filelist.txt: $!\n";
    foreach $sample (sort keys %{$samples}) {
	if (-e "$args->{scratchfolder}/singlesamples/$sample/fastqc/R1" . $ss . "_fastqc/fastqc_data.txt") {
	    print OFILE "$args->{scratchfolder}/singlesamples/$sample/fastqc/R1" . $ss . "_fastqc/fastqc_data.txt\t$sample\tR1\n";
	}
	if (-e "$args->{scratchfolder}/singlesamples/$sample/fastqc/R2" . $ss . "_fastqc/fastqc_data.txt") {
	    print OFILE "$args->{scratchfolder}/singlesamples/$sample/fastqc/R2" . $ss . "_fastqc/fastqc_data.txt\t$sample\tR2\n";
	}
    }
    close OFILE;
    $result = `fastqqualityplot.pl -q -d -l fastqqualityplot-filelist.txt`;
    `mv fastqqualityplot.png $args->{outputfolder}/fastqqualityplot.png`;
    print $result;

    ### Filtered fastqqualityplot ###
    if ($args->{qualitycontrol}) {
	print "Running filtered fastqqualityplot.pl\n" if ($args->{verbose});
	open OFILE, ">filtered_fastqqualityplot-filelist.txt" or die "cannot open filtered_fastqqualityplot-filelist.txt: $!\n";
	foreach $sample (sort keys %{$samples}) {
	    if (-e "$args->{scratchfolder}/singlesamples/$sample/fastqc-trim/R1" . $ss . ".trim_fastqc/fastqc_data.txt") {
		print OFILE "$args->{scratchfolder}/singlesamples/$sample/fastqc-trim/R1" . $ss . ".trim_fastqc/fastqc_data.txt\t$sample\tR1\n";
	    }
	    if (-e "$args->{scratchfolder}/singlesamples/$sample/fastqc-trim/R2" . $ss . ".trim_fastqc/fastqc_data.txt") {
		print OFILE "$args->{scratchfolder}/singlesamples/$sample/fastqc-trim/R2" . $ss . ".trim_fastqc/fastqc_data.txt\t$sample\tR2\n";
	    }
	}
	close OFILE;
	$result = `fastqqualityplot.pl -q -d -l filtered_fastqqualityplot-filelist.txt`;
	`mv fastqqualityplot.png $args->{outputfolder}/filtered_fastqqualityplot.png`;
	print $result;
    }

    ### dimer and adapter plot ###
    # generate reads per sample plot
    $ofile = "fastqcdimer.tmp";
    $rfile = "fastqcdimer.r";
    print "Generating dimer/adapter plot\n" if ($args->{verbose});
    open OFILE, ">$ofile" or die "cannot open temporary file $ofile for writing: $!\n"; 
    # print out the data
    print OFILE "sample\tdata\tvalue\n";
    foreach $sample (sort keys %{$samples}) {
	$dimer = $stats->{$sample}{fastqc}{"\%dimer"} // 0;
	$adapter = $stats->{$sample}{fastqc}{"\%adapter"} - $dimer;
	$good = 100 - $dimer - $adapter;
	print OFILE "$sample\tdimer\t$dimer\n";
	print OFILE "$sample\tadapter\t$adapter\n";
	print OFILE "$sample\tgood\t$good\n";
    }
    close OFILE;

    my $height = 480;
    my $width = 480;
    my $numsamples = keys %{$samples};
    if ($numsamples > 6) {
	$width = 480 + (20 * ($numsamples-6));
    }

    open RFILE, ">$rfile" or die "Cannot open $rfile\n";
    print RFILE qq(
  library(ggplot2);
  library(scales);
  datat <- read.table("$ofile", header=T, colClasses=c("sample"="factor"));
  png(filename="dimer.png", height = $height, width = $width);

  ggplot(datat, aes(x = sample, y = value, fill = data)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) + xlab("Sample") + ylab("Percent of Reads") + theme(legend.title=element_blank()) + scale_y_continuous(labels = comma) + scale_fill_manual(values=c("#B50000", "#7F7F7F", "#00008B"));

  dev.off();

#  png(filename="dimer.v2.png", height = $height, width = $width);

#  ggplot(datat, aes(reorder(sample, -value), value, fill=data)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) + xlab("Sample (sorted by Y-value)") + ylab("Number of Reads") + theme(legend.title=element_blank()) + scale_y_continuous(labels = comma) +scale_fill_manual(values=c("#7F7F7F", "#00008B"))

#  dev.off();

  #eof
);

    close RFILE;

    system("R --no-restore --no-save --no-readline < $rfile > $rfile.out");

    `mv *.png $args->{outputfolder}`;

    # Add to report
    &h1("Fastq Quality Plots");
    
    $collapse = (keys %{$samples} > 6) ? 1 : 0;
    $space = "";
    if ($collapse) {
	&line(".. container:: toggle");
	&line("  .. container:: header");
        &line("    **Show/Hide FastQC file list**");
	&line("  .. container:: FastQC Files");
	$space = "    ";
    }

    foreach $sample (sort keys %{$samples}) {
#	&bullet("$sample");
	print $reportFH "${space}- $sample `R1-FastQC <Resources/fastqc/$sample" . "_R1.html>`__";
	print $reportFH " `R2-FastQC <Resources/fastqc/$sample" . "_R2.html>`__" if ($args->{pe});
	print $reportFH " `Run-log <Resources/logs/$sample" . "-log.txt>`__";
	print $reportFH "\n";
    }
    &space;
    &image("fastqqualityplot.png");
    &image("dimer.png");
    &space;

    if ($args->{qualitycontrol}) {
	&h1("Filtered Fastq Quality Plots");

	if ($collapse) {
	    &line(".. container:: toggle");
	    &line("  .. container:: header");
	    &line("    **Show/Hide Filtered FastQC file list**");
	    &line("  .. container:: Filtered FastQC Files");
	    $space = "    ";
	}

	foreach $sample (sort keys %{$samples}) {
#	&bullet("$sample");
	    print $reportFH "${space}- $sample `R1-FastQC <Resources/fastqc_filtered/$sample" . "_R1.html>`__.";
	    print $reportFH " `R2-FastQC <Resources/fastqc_filtered/$sample" . "_R2.html>`__" if ($args->{pe});
	    print $reportFH "\n";
	}
	&space;
	&image("filtered_fastqqualityplot.png");
	&space;
    }

}

############################## Expressiontable Plot #######################
sub expressiontableplot {
    my ($samples, $args, $stats, $expressiontable, $name) = @_;

    print "\nGenerating Expressiontable Plot\n" if ($args->{verbose});

    ### expressiontableplot ###
    print "Running expressiontableplot.pl\n" if ($args->{verbose});
    @result = `expressiontableplot.pl $expressiontable`;
    `rm dendrogramplot.png`;
    `rm mdsplot.png`;
    # move and rename plots
    `for f in *.png ; do mv "\$f" "$args->{outputfolder}/${name}_\$f" ; done`;
    # save stats from stdout
    $statzone = 0;
    foreach $result (@result) {
#	print $result if ($args->{verbose});
	$statzone = 0 if ($result =~ /^STATS End/);
	if ($statzone) {
	    chomp $result;
	    ($sample, $key, $stat) = split /\t/, $result;
	    $stats->{$sample}{$name}{$key} = $stat;
	}
	$statzone = 1 if ($result =~ /^STATS Start/);
    }

    ### expressiontablepca ###
    print "Running expressiontablepca.pl\n" if ($args->{verbose});
    if ($args->{singlecell}) {
	$result = `expressiontablepca-test.pl $expressiontable $args->{samplesheet}`;
    } else {
	$result = `expressiontablepca.pl $expressiontable $args->{samplesheet}`;
    }
    # move and rename plots
    `for f in *.png ; do mv "\$f" "$args->{outputfolder}/${name}_\$f" ; done`;
    `mkdir $args->{outputfolder}/${name}_emperor`;
    `cp -r emperor $args->{outputfolder}/${name}_emperor`;
    `cp -r DESeq2-data.RDATA $args->{outputfolder}/${name}_DESeq2-data.RDATA`;

    # add to report
    &h2("$expressiontable plots");
    &line("The dendrogram, mdsplot, and interactive PCA plot are all generated from the top 500 most variable genes. The heatmap is generated from the top 50 most variable genes.");
    &image("${name}_expressedplot.png");
    &image("${name}_expressiondistributionplot.png");
    &image("${name}_dendrogramplot.png");
    &image("${name}_mdsplot.png");
    &image("${name}_heatmapplot.png");
    &line("|pcaicon| `Interactive PCA plot <Resources/${name}_emperor/emperor/index.html>`_"); 
    &line(".. |pcaicon| image:: Resources/pcaicon.png"); 
    &space;
}

############################## Metrics #######################
sub metrics {
    print "\nRunning Metrics\n";

    my ($samples, $args, $stats, $metricscolumns) = @_;
    my @columns = @{$metricscolumns};

    my $ofile = "$args->{outputfolder}/metrics.txt";
    open OFILE, ">$ofile" or print "Cannot open file $ofile for writing: $!\n";
    foreach $column (@columns) {
        print OFILE "$column\t";
        print "$column\t";
    }
    print "\n";
    print OFILE "\n";

    foreach $sample (sort keys %{$stats}) {
        foreach $column (@columns) {
            if (defined ($stats->{$sample}{$column})) {
                print $stats->{$sample}{$column} . "\t";
                print OFILE $stats->{$sample}{$column} . "\t";
            } else {
                print "undef\t";
                print OFILE "undef\t";
            }
        }
        print "\n";
        print OFILE "\n";
    }

    # calcluate mean values
    print "$args->{runname}    \tMean    \t";
    print OFILE "$args->{runname}    \tMean    \t";
    foreach $column (@columns[2..$#columns]) {
        $sum = 0;
        $count = 0;
        foreach $sample (sort keys %{$stats}) {
            if (defined ($stats->{$sample}{$column}) and looks_like_number($stats->{$sample}{$column})) {
                $sum += $stats->{$sample}{$column};
                $count++;
            }
        }
        if ($count > 0) {
            print round10($sum / $count);
            print "\t";
            print OFILE round10($sum / $count);
            print OFILE "\t";
        } else {
            print "undef\t";
            print OFILE "undef\t";
        }
    }
    print "\n";
    print OFILE"\n";

}

############################## allmetrics #######################
sub allmetrics {
    print "\nRunning All Metrics\n";

    my ($samples, $args, $stats) = @_;

    # open output file
    my $ofile = "$args->{outputfolder}/metrics.txt";
    open OFILE, ">$ofile" or print "Cannot open file $ofile for writing: $!\n";

    # get list of all stages and all stats of each stage
    foreach $sample (keys %{$samples}) {
	foreach $stage (keys %{$stats->{$sample}}) {
	    foreach $stat (keys %{$stats->{$sample}{$stage}}) {
		$stages{$stage}{$stat}++;
		$stagelist{$stage} = 1;
	    }
	}
    }
    
    # determine order for printing stages
    # use order recorded in $args->{stageorder}
    foreach $stage (@{$args->{stageorder}}) {
	if (defined($stagelist{$stage})) {
	    push @stageorder, $stage;
	    delete $stagelist{$stage};
	}
    }
    # add any additional stages on to the end
    foreach $stage (sort keys %stagelist) {
	next if ($stage eq "temp");
	next if ($stage eq "files");
	push @stageorder, $stage;
    }
    
    # print header
#    foreach $stage (sort keys %stages) {
    foreach $stage (@stageorder) {
	foreach $stat (sort keys %{$stages{$stage}}) {
	    print OFILE "$stage-$stat\t";
	    print "$stage-$stat\t";
	}
    }
    print "\n";
    print OFILE "\n";

    # print data
    foreach $sample (sort keys %{$stats}) {
#	foreach $stage (sort keys %stages) {
	foreach $stage (@stageorder) {
	    foreach $stat (sort keys %{$stages{$stage}}) {

		if (defined ($stats->{$sample}{$stage}{$stat})) {
		    chomp $stats->{$sample}{$stage}{$stat}; # need this chomp, don't know why 
		    print $stats->{$sample}{$stage}{$stat} . "\t";
		    print OFILE $stats->{$sample}{$stage}{$stat} . "\t";
		} else {
		    print "undef\t";
		    print OFILE "undef\t";
		}
	    }
	}
        print "\n";
        print OFILE "\n";
    }

    # calcluate mean values
    foreach $stage (@stageorder) {
#    foreach $stage (sort keys %stages) {
	foreach $stat (sort keys %{$stages{$stage}}) {

	    $sum = 0;
	    $count = 0;
	    foreach $sample (sort keys %{$stats}) {
		if (defined ($stats->{$sample}{$stage}{$stat}) and looks_like_number($stats->{$sample}{$stage}{$stat})) {
		    $sum += $stats->{$sample}{$stage}{$stat};
		    $count++;
		}
	    }
	    if (($stage eq "general") and ($stat eq "runname")) {
		print "$args->{runname}\t";
		print OFILE "$args->{runname}\t";
	    } elsif (($stage eq "general") and ($stat eq "samplename")) {
		print "Mean\t";
		print OFILE "Mean\t";
	    } elsif ($count > 0) {
		print round100($sum / $count);
		print "\t";
		print OFILE round100($sum / $count);
		print OFILE "\t";
	    } else {
		print "undef\t";
		print OFILE "undef\t";
	    }
	}
    }
    print "\n";
    print OFILE"\n";
    close OFILE;

    print "Scratch folder: $args->{scratchfolder}\n";
    print "Output folder: $args->{outputfolder}\n";

}

############################## Sphinx #######################
sub sphinx {
    print "\nRunning Sphinx to generate HTML report\n";

    my ($samples, $args, $stats, $reporttype, $reportfilename) = @_;

    my $outputfolder = $args->{outputfolder};
    my $scratchfolder = $args->{sphinxscratchfolder};
    my $myscratchfolder = "$args->{scratchfolder}/sphinx";
    my $log = "$args->{scratchfolder}/sphinxlog";
    my $htmlfile = "$outputfolder/index.html";

    # copy over sphinx template copy over rst files (metrics, analysisreport, DESeq2, etc)
    `cp -r $ENV{GOPHER_PIPELINES}/resources/sphinx_directory $myscratchfolder`;
    `sed -i "s|REPORTTITLE|$reporttype|" $myscratchfolder/conf.py`;

    # symlink outputfolder contents over
    `ln -s $outputfolder/* $myscratchfolder`;

    # generate html report with Sphinx
    chdir $myscratchfolder;
    $result = `module load python/2.7.1; make html &> $log`;
    print $result if ($args->{verbose});
    $result = `cat $log`;
    print $result if ($args->{verbose});
    `mv $log $myscratchfolder`;

    # Move outputfolder to Resources folder, copy sphinx output over
    `mkdir $outputfolder/Resources`;
    `mv $outputfolder/* $outputfolder/Resources &> /dev/null`;
    `cp -r $myscratchfolder/_build/html/* $outputfolder/Resources`;

    # rename html file and update links
    `cp $outputfolder/Resources/AnalysisReport.html $htmlfile`;
    # the _images folder is buried in the Resources folder - adjust links
    `sed -i "s|_images/|Resources/|g" $htmlfile`;
    # the _static folder is buried in the Resources folder - adjust links
    `sed -i "s|_static/|Resources/_static/|g" $htmlfile`;
    # enable javascript to swap images on click
    `sed -i "s|<img |<img onclick=imageswap(this) |g" $htmlfile`;
    `sed -i "s|</head>|</title>\\n<script>\\nfunction imageswap(x) { if (x.src.endsWith('v1.png')) { x.src = x.src.replace('v1.png','v2.png'); } else if (x.src.endsWith('v2.png')) { x.src = x.src.replace('v2.png', 'v1.png'); } else { x.src = x.src; } }\\n</script>\\n</head>|" $htmlfile`;

    # copy over programversions and metrics - this is a hack
#    `cp $myscratchfolder/*.txt $outputfolder/Resources`;

}

############################## Sphinx #######################
# Requires: sphinx-quickstart
sub oldsphinx {
    print "\nRunning Sphinx to generate HTML report\n";

    my ($samples, $args, $stats, $reporttype, $reportfilename) = @_;

    my $outputfolder = $args->{outputfolder};
    my $scratchfolder = $args->{scratchfolder};
    my $myscratchfolder = "$args->{scratchfolder}/sphinx";
    my $log = "$args->{scratchfolder}/sphinxlog";
    my $htmlfile = "$outputfolder/../$reportfilename";

    # generate html report with Sphinx
    `rm -r $myscratchfolder` if (-e $myscratchfolder); # sphinx wants to create the output folder
    # create sphinx directory, copy over rst files (metrics, analysisreport, DESeq2, etc)
    $result = `module load python/2.7.1; sphinx-quickstart --quiet --project="$reporttype" --author="Regents of the University of Minnesota" --master="AnalysisReport" -v=1 $myscratchfolder > $log; cd $myscratchfolder; cp $outputfolder/*.rst $myscratchfolder; cp $outputfolder/metrics.txt $myscratchfolder; cp $outputfolder/programversions.txt $myscratchfolder`;
    chdir $myscratchfolder;
    mkdir "_themes";
    `ln -s $ENV{GOPHER_PIPELINES}/resources/sphinx_rtd_theme $myscratchfolder/_themes/sphinx_rtd_theme`;

    # put images in a local directory to avoid warnings from sphinx
    `mkdir $myscratchfolder/Resources`;
    `ln -s $outputfolder/*.png $myscratchfolder/Resources/`;
    `ln -s $outputfolder/*.gif $myscratchfolder/Resources`;
    `ln -s $outputfolder/*.jpg $myscratchfolder/Resources`;

    `cp $ENV{GOPHER_PIPELINES}/resources/sphinx_files/page.html $myscratchfolder/_templates`;
    `cp $ENV{GOPHER_PIPELINES}/resources/sphinx_files/custom.css $myscratchfolder/_static`;
    
    `sed -i "s|html_theme = 'default'|html_theme = 'sphinx_rtd_theme'|g" $myscratchfolder/conf.py`;
    $command = 'sed -i "s|#html_theme_path = \[\]|html_theme_path = \[\'_themes\',\]|g" ' . "$myscratchfolder/conf.py";
    `$command`;

    $result = `module load python/2.7.1; make html &> $log`;
    print $result if ($args->{verbose});
    $result = `cat $log`;
    print $result if ($args->{verbose});
    `mv $log $myscratchfolder`;
    # move html file to output folder
    `cp $myscratchfolder/_build/html/AnalysisReport.html $htmlfile`;
    # the _images folder is buried in the Resources folder - adjust links
    `sed -i "s|_images/|Resources/|g" $htmlfile`;
    # the _static folder is buried in the Resources folder - adjust links
    `sed -i "s|_static/|Resources/_static/|g" $htmlfile`;
    # enable javascript to swap images on click
    `sed -i "s|<img |<img onclick=imageswap(this) |g" $htmlfile`;
    `sed -i "s|</head>|</title>\\n<script>\\nfunction imageswap(x) { if (x.src.endsWith('v1.png')) { x.src = x.src.replace('v1.png','v2.png'); } else if (x.src.endsWith('v2.png')) { x.src = x.src.replace('v2.png', 'v1.png'); } else { x.src = x.src; } }\\n</script>\\n</head>|" $htmlfile`;
    # move the rest of the html supporting files
    `cp -r $myscratchfolder/_build/html/* $outputfolder`;
}

# make sure columns are tab delimited - convert multiple spaces and tabs to a single tab, allow the last column to contain spaces
sub fixsamplesheet {
    my ($samplesheet) = @_;
 
    # convert mac and dos newlines to unix newlines
    `perl -i -pe 's/\015/\012/g' $samplesheet`;
    `perl -i -pe 's/\015\012/\012/g' $samplesheet`;
   
    open IFILE, "$samplesheet" or die "Cannot open samplesheet $samplesheet: $!\n";
    open OFILE, ">samplesheet.tmp" or die "Cannot open temporary file samplesheet.tmp: $!\n";
    $header = <IFILE>;
    $header =~ s/^\s+|\s+$//g; # remove leading and trailing whitespace
    @header = split ' ', $header;
    $columns = $#header + 1;
    print OFILE join("\t", @header);
    print OFILE "\n";
    while ($line = <IFILE>) {
	$line =~ s/^\s+|\s+$//g; # remove leading and trailing whitespace
	next if ($line eq ""); # skip blank lines
	$line =~ s/ +\t/\t/g; # collapse spaces and tab into a tab
	$line =~ s/\t +/\t/g; # collapse tab and spaces into a tab
	@line = split /[ +|\t]/, $line, $columns;
	print OFILE join("\t", @line);
	print OFILE "\n";
    }
    $result = `mv samplesheet.tmp $samplesheet`;
    print $result;
    close OFILE;
    close IFILE;

}

############################################################################
################################ Helper subs ###############################
############################################################################

### die gracefully ###
sub diegracefully {
# catch interrupt signals
    $SIG{'INT'} = sub { &runtime; print "WARNING: program exited prematurely, results are incomplete\n"; exit; };
    $SIG{'TERM'} = sub { &runtime; print "WARNING: program exited prematurely, results are incomplete\n"; exit; };
    $SIG{'KILL'} = sub { &runtime; print "WARNING: program exited prematurely, results are incomplete\n"; exit; };
}

### Read in extra options file ###
sub readinextraoptions {
    my ($args) = @_;

    if ($args->{extraoptionsfile}) {
        open $fh, $args->{extraoptionsfile} or die "Cannot open extra command file $args->{extraoptionsfile}: $!\n";
        while ($line = <$fh>) {
            chomp $line;
            next if ($line =~ /^#/);
            next if ($line eq "");
            ($command, $options) = split /\s+/, $line, 2;
            $command = lc($command);
            print "Using extra options for command $command: $options\n";
            $args->{extraoptions}{$command} = $options;
        }
    }
}

# Round to neaerest tenth
sub round10 {
    return sprintf("%.1f", $_[0]);
#    return int($_[0] * 10 + 0.5) / 10;
}

# Round to neaerest hundreth
sub round100 {
    return sprintf("%.2f", $_[0]);
#    return int($_[0] * 100 + 0.5) / 100;
}

sub printlog {
    my ($text, $logfile) = @_;
    open LOG, ">>$logfile" or die "Cannot open log file $logfile: $!\n";
    print LOG "# $text\n";
    print "$text\n";
}

# run a system command, pipe stderr to stdout, append stdout to a log file (if defined), and return stdout as one big string
sub run {
    my ($command, $logfile, $output) = @_;

    my ($result, $line);

#    my $redirect = "2>\%1";
#    if ($output) {
#	if ($output eq "stderr") {
#	    $redirect = 

    $loglines = 0;
    if (defined($logfile)) {
	open LOG, ">>$logfile" or die "Cannot open log file $logfile: $!\n";
	print LOG "\n#################################################################\n";
	print LOG "COMMAND: $command\n";
	open COMMAND, "$command 2>\&1 |" or warn "Cannot execute $command\n$!\n";
	$result = "";
	while ($line = <COMMAND>) {
	    $result .= $line;
	    $loglines++;
	    print LOG "Output truncated at 1000 lines\n" if ($loglines eq 1000);
	    next if ($loglines >= 1000);
	    print LOG $line;

	}
	close COMMAND;
    } else {
	print "\n#################################################################\n";
	print "COMMAND: $command\n";
	$result = `$command 2>\&1`;
	print $result;
    }
    return $result;
}

# print how long we've been running
sub runtime {
    state $firsttime = 1;
    state $starttime = time();
    state $checkpointtime = $starttime;
    $now = time();
    my $runtime = $now - $starttime;
    my $checktime = $now - $checkpointtime;
    if ($firsttime) {
	$checkpointtime = $now;
	$firsttime = 0;
	return;
    }
    print "Checkpoint: " . &formattime($checktime) . "\tTotal time: " . &formattime($runtime) . "\n";
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

# determine if a file has a .gz extension
sub compressed {
    if ($_[0] =~ /\.gz$/) {
	return 1;
    } else {
	return 0
    }
}

# function to run when an interrupt signal has been received
sub INT_handler {

    # print out a report
    &runtime;

    print "WARNING: program exited prematurely, results are incomplete\n";

    exit;
}

# return the memory available on this node in GB
sub nodememory {

    my $info = `head -n 1 /proc/meminfo | head -n 1`;
    my ($junk, $mem_info) = split /\s+/, $info;
    chomp $mem_info;
    $mem_info = $mem_info / 1024; # turn into mb
    $mem_info = int($mem_info / 1024 * 10 + .5) / 10; # turn into gb
    return $mem_info;

}

### check that required software is present
sub requiredprograms {
    my @programs = @_;

    my $quit = 0;
    foreach $program (@programs) {
        # check if shell variables are defined (good for checking for java software)
#	print "Checking for $program\n";
        if ($program =~ /^\$/) {
            $result = `echo $program`;
            chomp $result;
            if ($result eq "") {
                print "Shell variable $program must be defined, please load the neccessary module or install the software\n";
                $quit = 1;
            }
	} else {
	    $result = `which $program 2>&1`;
#	    print "$result\n";
	    if ($result =~ "no $program in") {
		print "Required program \"$program\" not found, please load the neccessary module or install the software\n";
		$quit = 1;
	    }
	}
    }
    exit if ($quit);
}

### Report versions of software used in the pipeline
sub programversions {
    my ($args, $programs) = @_;

    # report gopher-pipelines version
    $text = "### Gopher-Pipelines Version ###";
    $command = "changelog | head"; 
    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";

    foreach $program (@{$programs}) {
	if (lc($program) eq "fastqc") {
	    $text = "### FastQC Version ###";
	    $command = "fastqc -v"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "fastq_screen") {
	    $text = "### Fastq_screen Version ###";
	    $command = "fastq_screen -version"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "prinseq") {
	    $text = "### Prinseq Version ###";
	    $command = "prinseq-lite.pl --version"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "trimmomatic") {
	    $text = "### Trimmomatic Version ###";
	    $command = "echo \$TRIMMOMATIC"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "samtools") {
	    $text = "### Samtools Version ###";
	    $command = "samtools --version"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "picard") {
	    $text = "### Picard Tools Version ###";
	    $command = "java -jar \$PICARDTOOLS -version 2>&1"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "cufflinks") {
	    $text = "### Cufflinks Version ###";
	    $command = "cufflinks 2>&1";
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "hisat2") {
	    $text = "### Hisat2 Version ###";
	    $command = "hisat2 --version"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "tophat2") {
	    $text = "### Tophat2 Version ###";
	    $command = "tophat2 --version"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "bowtie2") {
	    $text = "### Bowtie2 Version ###";
	    $command = "bowtie2 --version"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "featurecounts") {
	    $text = "### Subread featureCounts Version ###";
	    $command = "featureCounts -v"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "bwa") {
	    $text = "### BWA Version ###";
	    $command = "bwa 2>&1";
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "pandaseq") {
	    $text = "### PandaSeq Version ###";
	    $command = "pandaseq -v 2>&1"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "cutadapt") {
	    $text = "### Cutadapt Version ###";
	    $command = "cutadapt --version"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "samstat") {
	    $text = "### Samstat Version ###";
	    $command = "samstat"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "salmon") {
	    $text = "### Salmon Version ###";
	    $command = "salmon --version 2>&1"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "kallisto") {
	    $text = "### Kallisto Version ###";
	    $command = "kallisto"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "r") {
	    $text = "### R Version ###";
	    $command = "R --version"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "qiime") {
	    $text = "### Qiime Version ###";
	    $command = "print_qiime_config.py"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "kraken") {
	    $text = "### Kraken Version ###";
	    $command = "kraken --version"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "krona") {
	    $text = "### Krona Version ###";
	    $command = "which ktextractTaxonomy"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "sstacks") {
	    $text = "### Stacks Version ###";
	    $command = "which $program"; 
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} elsif (lc($program) eq "tassel") {
	    $text = "### Tassel Version ###";
	    $command = "run_pipeline.pl | grep Version";
	    $return .= $text . "\n" . $command . "\n" . `$command` . "\n\n";
	} else {
	    $text = "### $program Version ###";
	    $command = "ERROR: this pipeline doesn't know how to identify the version of this program\n";
	    $return .= $text . "\n" . $command . "\n" . `which $program` . "\n\n";
	}
    }
    $ofile = "$args->{outputfolder}/programversions.txt";
    open OFILE, ">$ofile" or die "Cannot open program version file $ofile: $!\n";
    print OFILE $return;
    close OFILE;
    return;

}

### Report generation subs ###

sub title {
    $text = shift @_;
    print $reportFH '=' x length($text) . "\n";
    print $reportFH "$text\n";
    print $reportFH '=' x length($text) . "\n";
}

sub subtitle {
    &h($_[0], "=");
}

sub h1 {
    &h($_[0], "-");
}

sub h2 {
    &h($_[0], ".");
}

sub h3 {
    &h($_[0], "_");
}

sub h {
    $text = shift @_;
    $char = shift @_;
    print $reportFH "$text\n";
    print $reportFH $char x length($text) . "\n";
}

sub line {
    $text = shift @_;
    print $reportFH "$text\n\n";
}

sub fieldlist {
    $text1 = shift @_;
    $text2 = shift @_;
    print $reportFH ":$text1: $text2\n";
}

sub bullet {
    $text = shift @_;
    print $reportFH "- $text\n";
}

sub space {
    print $reportFH "\n";
}

sub image {
    $file = shift @_;
    print $reportFH ".. image:: $file\n";
}

# inline link
sub link {
    $url = shift @_;
    $text = shift @_;
    print $reportFH "`$text <$url>`_";
}


############ Debugging helper subs ##########
sub statprint {
    my ($samples, $args, $stats) = @_;
    print "DEBUG STATS:\n";
    foreach $sample (sort keys %{$stats}) {
	print "$sample";
	foreach $column (sort keys %{$stats->{$sample}}) {
	    print "\t$column";
	}
	print "\n";
    }
    print "\n";
}


# end of Pipline Module
1;
