#!/usr/bin/perl -w

######################################################
# rnaseq-pipeline.pl
# John Garbe
# December 2015
#
# TODO: add option for controlling spike bowtie
#
# Get rid of extra blank pages in pdf sphinx output:
# Put this in your source/conf.py in the "Options for LaTeX output" section:
# 'classoptions': ',openany,oneside', 'babel' : '\\usepackage[polish]{babel}'
#
#######################################################

=head1 DESCRIPTION

rnaseq-pipeline.pl - RNA-seq Analysis Program: Analyze a set of RNA-seq fastq files using Hisat2 (or Tophat2) and Cufflinks

=head1 SYNOPSIS

rnaseq-pipeline.pl --hisat2index index --gtffile file --fastqfolder folder
rnaseq-pipeline.pl --bowtie2index index --gtffile file --fastqfolder folder

=head1 OPTIONS

Analyze a collection of rna-seq fastq files

 --outputfolder folder : A folder to deposit final results
 --samplespernode integer : Number of samples to process simultaneously on each node (default = 1)
 --threadspersample integer : Number of threads used by each sample
 --subsample integer : Subsample the specified number of reads from each sample. 0 = no subsampling (default = 0)
 --qualitycontrol : Enable quality-control: use trimmomatic to trim adapter sequences, trim low quality bases from the ends of reads, and remove short sequences
 --extraoptionsfile file : File with extra options for trimmomatic, tophat, cuffquant, or featurecounts
 --fastqfolder folder : A folder containing fastq files to process
 --hisat2index index : A hisat2 index basename (for running hisat2)
 --bowtie2index index : A bowtie2 index basename (for running tophat2)
 --gtffile file : A reference genome gtf gene annotation file
 --maskfile : GTF file for Cuffquant's --maskfile option
 --stranded : TruSeq Stranded dataset
 --samplesheet file : A samplesheet
 --name string : Name of the analysis, fastqfolder name is default name
 --scratchfolder folder : A temporary/scratch folder
 --singlecell : Single-cell analysis (uses a different data normalization before generating PCA plots)
 --spikefasta file : A spike-in control fasta file
 --headcrop integer : crop the specified number of bases off the front of each read, useful for picogreen libraries
 --resume : Continue where a failed/interrupted run left off
 --help : Print usage instructions and exit
 --verbose : Print more information while running

=cut

##################### Main ###############################

use FindBin;                  # locate this script
use lib "$FindBin::RealBin";  # use the current directory
use Pipelinev4;

use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use Pod::Usage;

my %samples;
my %args;
my %stats;
my $reportFH; # analysis report filehandle

### Initialize ###
&init(\%samples, \%args, \%stats);

### Start report ###
startreport(\%samples, \%args, \%stats, "RNA-Seq Analysis Report");

### Singlesample analysis ###
&singlesamples(\%samples, \%args, \%stats);
runtime;

### Cuffnorm ###
&cuffnorm(\%samples, \%args, \%stats);
runtime;

### Differential expression ###
#&differentialexpression(\%samples, \%args, \%stats);
runtime;

### All sample processing ###
&allsample(\%samples, \%args, \%stats);
runtime;

### Finish report ###
&finishreport(\%samples, \%args, \%stats);
runtime;

### Metrics ###
runtime;
allmetrics(\%samples, \%args, \%stats);

### Sphinx ###
sphinx(\%samples, \%args, \%stats, "RNA-Seq Analysis", "index.html");
runtime;

print "Finished\n";
runtime;

exit;

##################### Initialize ###############################
sub init {
    print "\nRunning Initialize\n";

    my ($samples, $args, $stats) = @_;

    # set defaults
    $args->{options} = join " ", @ARGV;
    my $help = "";
    $args->{threadspersample} = 8;
    if ($ENV{"PBS_NUM_PPN"}) {
	$args->{samplespernode} = int($ENV{"PBS_NUM_PPN"} / $args->{threadspersample} + 0.5);
    } else {
	$args->{samplespernode} = 1;
    }
    GetOptions("help" => \$args->{help},
               "verbose" => \$args->{verbose},
               "resume" => \$args->{resume},
               "threadspersample=i" => \$args->{threadspersample},
               "samplespernode=i" => \$args->{samplespernode},
               "outputfolder=s" => \$args->{outputfolder},
               "scratchfolder=s" => \$args->{scratchfolder},
               "extraoptionsfile=s" => \$args->{extraoptionsfile},
               "subsample=i" => \$args->{subsample},
               "fastqfolder=s" => \$args->{fastqfolder},
               "samplesheet=s" => \$args->{samplesheet},
               "name=s" => \$args->{runname},
               "qualitycontrol" => \$args->{qualitycontrol},
               "gtffile=s" => \$args->{gtffile},
               "maskfile=s" => \$args->{maskfile},
               "hisat2index=s" => \$args->{hisat2index},
               "bowtie2index=s" => \$args->{bowtie2index},
               "singlecell" => \$args->{singlecell},
               "spikefasta=s" => \$args->{spikefasta},
               "stranded" => \$args->{stranded},	       
	       "headcrop=i" => \$args->{headcrop},
        ) or pod2usage;
    pod2usage(-verbose => 99, -sections => [qw(DESCRIPTION|SYNOPSIS|OPTIONS)]) if ($args->{help});
    pod2usage unless ( ($args->{hisat2index} or $args->{bowtie2index})  and $args->{fastqfolder} and $args->{gtffile} );
    if ($#ARGV >= 0) {
	print "Unknown commandline parameters: @ARGV\n";
	pod2usage;
    }
    $args->{threads} = $ENV{"PBS_NUM_PPN"} // $args->{threadspersample}; # threads used by pipeline stages (not singlesample stages)
    $args->{extraoptionsfile} = abs_path($args->{extraoptionsfile}) if ($args->{extraoptionsfile});

    ### Starttime ###
    my $date = `date`;
    chomp $date;
    print "Starting at $date\n";
    $args->{startdate} = $date;
    &runtime;

    ### Handle Parameters ###
    # add parameter processing here

    # $args->{spikefasta} = "/panfs/roc/itascasoft/rnaseqQC/1.1/data/spike.fasta";
    die "ERROR: Could not locate gtffile $args->{gtffile}\n" if (! -e $args->{gtffile});
    $args->{gtffile} = abs_path($args->{gtffile});
    if ($args->{maskfile}) {
	die "ERROR: Could not locate maskfile $args->{maskfile}\n" if (! -e $args->{maskfile});
	$args->{maskfile} = abs_path($args->{maskfile});
    }
    if ($args->{hisat2index}) {
	$args->{hisat2index} = hisat2indexcheck($args->{hisat2index});
    } else {
	$args->{bowtie2index} = bowtie2indexcheck($args->{bowtie2index});
    }

    ### Finish setup
    requiredprograms(("fastqc", "\$TRIMMOMATIC", "\$PICARD", "parallel", "R", "samtools", "cufflinks", "sphinx-quickstart")); # check that required programs are installed
    diegracefully();
    readinextraoptions($args);
    samplesheetandfoldersetup($samples, $args, $stats, "rnaseq");
    my @programs = ("fastqc", "trimmomatic", "picard", "R", "samtools", "cufflinks"); # record versions of programs used
    if ($args->{hisat2index}) {
	push @programs, "hisat2";
    } else {
	push @programs, "tophat2";
    }
    programversions($args, \@programs);

}

######################## Singlesample jobs ##########################
sub singlesamples {
    print "\nRunning Singlesample jobs\n";

    my ($samples, $args, $stats) = @_;

    my $myscratchfolder = "$args->{scratchfolder}/singlesamples";
    mkdir $myscratchfolder;

    # set up nodefile
    $nodefile = "$myscratchfolder/nodelist.txt";
    $nodecount = 1;
    my $junk;
    if ($ENV{PBS_NODEFILE}) {
	$result = `sort -u $ENV{PBS_NODEFILE} > $nodefile`;
	print $result;
	$nodecount = `wc -l < $nodefile`;
	chomp $nodecount;
    }

    # pass arguments through to singlesample.pl
    my $qc = $args->{qualitycontrol} ? "--qualitycontrol" : "";
    my $subsample = $args->{subsample} ? "--subsample $args->{subsample}" : "";
    my $spike = $args->{spikefasta} ? "--spikefata $args->{spikefasta}" : "";
    my $extraoptions = $args->{extraoptionsfile} ? "--extraoptionsfile $args->{extraoptionsfile}" : "";
    my $stranded = $args->{stranded} ? "--stranded" : "";
    my $mask = $args->{maskfile} ? "--maskfile $args->{maskfile}" : "";
    my $headcrop = $args->{headcrop} ? "--headcrop $args->{headcrop}" : "";
    if ($args->{hisat2index}) {
	$index = "--hisat2index $args->{hisat2index}";
    } else {
	$index = "--bowtie2index $args->{bowtie2index}";
    }
    # print out singlesample commandlines
    my $commandfile = "$myscratchfolder/singlesample-commands.txt";
    open OFILE, ">$commandfile" or die "Cannot open $commandfile: $!\n";
    foreach $sample (keys %{$samples}) {
	if (($args->{resume}) && (-e "$myscratchfolder/$sample/Complete")) {
	    print "Skipping completed sample $sample\n";
	    next;
	}
	if ($args->{pe}) {
	    $command = "rnaseq-singlesample.pl --threads $args->{threadspersample} $qc $subsample $stranded $mask --outputfolder $myscratchfolder/$sample --gtffile $args->{gtffile} $index $spike $headcrop $extraoptions --R1 " . $samples->{$sample}{R1}{fastq} . " --R2 " . $samples->{$sample}{R2}{fastq} . "\n";
	} else {
	    $command = "rnaseq-singlesample.pl --threads $args->{threadspersample} $qc $subsample $stranded $mask --outputfolder $myscratchfolder/$sample --gtffile $args->{gtffile} $index $spike $headcrop $extraoptions --R1 " . $samples->{$sample}{R1}{fastq} . "\n";
	}
	print OFILE $command;
    }
    close OFILE;

    # run jobs in parallel
    if ($nodecount <= 1) { # single node
	run("cat $commandfile | parallel -j $args->{samplespernode}", $args->{logfile});
    } else { # multiple node
	run("cat $commandfile | parallel -j $args->{samplespernode} --sshloginfile $nodefile --workdir $ENV{PWD}", $args->{logfile});
    }

    # read in stats from each sample
    getsinglesamplestats(\%samples, \%args, \%stats);

}

######################## Cuffnorm ##########################
sub cuffnorm {
    print "\nRunning Cuffnorm\n";

    my ($samples, $args, $stats) = @_;

    my $myscratchfolder = "$args->{scratchfolder}/cuffnorm";
#    mkdir $myscratchfolder;

    # generate list of cxb files
    @samples = sort keys %{$samples};
    my $cxbfiles = "$args->{scratchfolder}/singlesamples/$samples[0]/abundances.cxb";
    for $i (1..$#samples) {
        $cxbfiles = "$cxbfiles $args->{scratchfolder}/singlesamples/$samples[$i]/abundances.cxb";
    }

    # run cuffnorm
    run("cuffnorm --no-update-check --output-dir $myscratchfolder --quiet -p $args->{threads} $args->{gtffile} $cxbfiles", $args->{logfile});

    # rename samplenames
    open OFILE, ">$myscratchfolder/cufftablerename-samplelist.txt" or die "cannot open $myscratchfolder/cufftablerename-samplelist.txt: $!\n";
    for $i (0..$#samples) {
	$number = $i+1;
        print OFILE "q$number\t$samples[$i]\n";
    }
    close OFILE;

    $result = `cufftablerename.pl $myscratchfolder/genes.fpkm_table $myscratchfolder/cufftablerename-samplelist.txt > $myscratchfolder/genes.fpkm_table.renamed`;
    print $result if ($args->{verbose});
#    `mv *.RDATA $outputfolder`;
    `cp $myscratchfolder/genes.fpkm_table.renamed $args->{outputfolder}`;
    `cp $myscratchfolder/isoforms.fpkm_table $args->{outputfolder}`;

}

################### Differential Gene Expression #######################
sub differentialexpression {
    print "Running differential expression\n";

    my ($samples, $args, $stats) = @_;

    my $myscratchfolder = "$args->{scratchfolder}/differentialexpression";
    mkdir $myscratchfolder;

    # get sample groups
    foreach $sample (keys %{$samples}) {
        if ($samples->{$sample}{group}) {
            push @groups, $samples->{$sample}{group};
            push @{$groups{$samples->{$sample}{group}}}, $sample;
        }
    }

    # quit if improper number of groups
    $groupnumber = keys %groups;
    if (($groupnumber < 2) or ($groupnumber > 3)) {
        if ($groupnumber == 0) {
            print "No groups defined in samplesheet, skipping differential expression testing\n";
        } else {
            print "$groupnumber groups defined in samplesheet, differential expression testing is only run if 2 or 3 groups are present\n";
        }
        return;
    }

    ### cuffdiff ###
    # generate command line
    $cuffdifffiles = "";
    foreach $group (keys %groups) {
        foreach $sample (@{$groups{$group}}) {
            $cuffdifffiles .= "$args->{scratchfolder}/singlesamples/$sample/abundances.cxb,";
        }
        chop $cuffdifffiles; # remove trailing ,
        $cuffdifffiles .= " ";
    }
    # run cuffdiff
    run("cuffdiff --output-dir $myscratchfolder --num-threads $args->{threads} $args->{gtffile} $cuffdifffiles", $args->{logfile});

    ### DESeq2 ###
    # TODO...

}

################### All-sample jobs #######################
sub allsample {
    print "\nRunning Allsample\n";

    my ($samples, $args, $stats) = @_;

    my $scratchfolder = $args->{scratchfolder};
    my $myscratchfolder = "$scratchfolder/allsamples";
    mkdir $myscratchfolder;
    my $outputfolder = $args->{outputfolder};
    chdir $myscratchfolder;

    @samples = sort keys %{$samples};

    ### generate plots ###
    subsampleplot($samples, $args, $stats);

    fastqcplot($samples, $args, $stats);

    if ($args->{hisat2index}) {
	hisat2plot($samples, $args, $stats);
    } else {
	tophat2plot($samples, $args, $stats);
    }	

    insertsizeplot($samples, $args, $stats);

    ### create folders ###
    # move singlesample logs to log folder
    compilefolder($samples, $args, $stats, "logs", "log", "-log.txt");
    # copy log files to the output folder
    `cp -rL $args->{scratchfolder}/logs $args->{outputfolder}`;

    # move accepted_hits.bam and .bam.bai to bam folder
    if ($args->{hisat2index}) {
	compilefolder($samples, $args, $stats, "bam", "hisat2.sort.bam", ".bam");
	compilefolder($samples, $args, $stats, "bam", "hisat2.sort.bam.bai", ".bam.bai");
    } else {
	compilefolder($samples, $args, $stats, "bam", "tophat_out/accepted_hits.sort.bam", ".bam");
	compilefolder($samples, $args, $stats, "bam", "tophat_out/accepted_hits.sort.bam.bai", ".bam.bai");
    }

    # cuffquant folder of .cxb files
    compilefolder($samples, $args, $stats, "cuffquant", "abundances.cxb", ".cxb");
    # save a tarball to outputfolder
    $result = `tar -zcf $args->{outputfolder}/cuffquant-cxb-files.tar.gz $args->{scratchfolder}/cuffquant/*.cxb`;
    print $result;

    ### subread merge ###
    print "Running subread-merge.pl\n" if ($args->{verbose});
    open OFILE, ">subread-samplelist.txt" or die "cannot open subread-samplelist.txt: $!\n";
    for $sample (@samples) {
        print OFILE "$scratchfolder/singlesamples/$sample/subread-counts.txt\t$sample\n";
    }
    close OFILE;
    $result = `subread-merge.pl -f subread-samplelist.txt > $myscratchfolder/subread.txt`;
    print $result;
    `cp $myscratchfolder/subread.txt $outputfolder`;

    ### Generate per-sample count files for Sue Rathe ###
    my $countdir = "$myscratchfolder/subread/";
    mkdir $countdir;
    $header = `head -n 1 $myscratchfolder/subread.txt`;
    chomp $header;
    @header = split /\t/, $header;
    for $i (1..$#header) {
	$column = $i + 1;
	`cut -f 1,$column $myscratchfolder/subread.txt > $countdir/$header[$i].txt`;
    }
    `tar -zcf $outputfolder/subread.tar.gz -C $myscratchfolder subread`;

    expressiontableplot(\%samples, \%args, \%stats, "subread.txt", "subread");
    `cp $outputfolder/genes.fpkm_table.renamed .`; # so only the filename shows up in the report, not the full path
    expressiontableplot(\%samples, \%args, \%stats, "genes.fpkm_table.renamed", "cuffquant");
    
}

############################ Finish Report #############################
sub finishreport {
    print "Finishing Report\n";

    my ($samples, $args, $stats) = @_;

    # copy additional files to output
    `cp $GPRESOURCES/metrics.rst $args->{outputfolder}`;
    `cp $GPRESOURCES/umnlogo.gif $args->{outputfolder}`;
    `cp $GPRESOURCES/UMII-full-banner.png $args->{outputfolder}`;
    `cp $GPRESOURCES/pcaicon.png $args->{outputfolder}`;
    `cp $GPRESOURCES/DESeq2.rst $args->{outputfolder}`;

    &h1("Data");
    &line("The primary output of this analysis pipeline is the quality-control information shown in this document, as well as expression data available in the following files.");
    &line("Gene expression tables, suitable for import into R packages such as DESeq, EdgeR, Singular (single-cell), and others for statistical analysis of gene expression. Most packages require raw count data and cannot use FPKM data:");
    &line("- Raw count expression table: `subread.txt <Resources/subread.txt>`_");
    &line("- FPKM expression table: `genes.fpkm_table.renamed <Resources/genes.fpkm_table.renamed>`_");
    &line("- Raw count expression data in DESeq2 format: `DESeq2-data.RDATA <Resources/subread_DESeq2-data.RDATA>`_  |  `Instructions for using this file <Resources/DESeq2.html>`_");
    &line("CXB expression files, suitable for use with Cuffdiff for statistical analysis of simple experimental designs:");
    &line("- CXB expression files: `cuffquant-cxb-files.tar.gz <Resources/cuffquant-cxb-files.tar.gz>`_");

    # print example cuffdiff command
    my $cuffdifffiles = "";
    foreach $sample (sort keys %{$samples}) {
        $cuffdifffiles .= "$sample.cxb,";
    }
    chop $cuffdifffiles; # remove last comma
    &line("To run Cuffdiff uncompress the cuffquant-cxb-files.tar.gz file and run the following command in the cuffquant folder. Cxb files in the same experimental group must be listed together and separated by a comma; a space must separate groups::\n\n\tcuffdiff --num-threads 1 $args->{gtffile} $cuffdifffiles\n");

    &h1("Intermediate and temporary files");
    my $deletedate = `date -d "30 days" +"%a %b %d %Y"`;
    chomp $deletedate;
    &line("While the primary output files contain the most important results of the analysis, some users may be interested in the intermediate files generated during the analysis of the experiment. Full output from this analysis is available for a limited time on the MSI filesystem at $args->{scratchfolder}. **This data will be deleted on $deletedate**. Copy the full output to your home directory::\n\n\tcp -r $args->{scratchfolder} ~/\n\nCopy an individual folder to your home directory::\n\n\tcp -rL $args->{scratchfolder}/fastq ~/\n");
    
    &line("Contents of $args->{scratchfolder}:");
    # all pipelines have these:
    &fieldlist("allsamples", "Temporary files used to create summary plots");
    &fieldlist("fastq", "Symlinks (shortcuts) to fastq files for all samples");
    &fieldlist("output", "Summary plots, html report, expression tables, fastqc output for all samples, and Cuffquant cxb files for all samples");
    &fieldlist("singlesamples", "Contains one folder for each sample containing the output from processing one sample");
    &fieldlist("logs", "log file for each sample");
    &fieldlist("sphinx", "Temporary files used to create html report");
    # unique to this pipeline:
    &fieldlist("bam", "Hisat2- (or Tophat2) generated bam alignment files for all samples");
    &fieldlist("fastqc", "FastQC output for all samples");
    &fieldlist("cuffnorm", "Full set of output files from Cuffnorm");
    &fieldlist("cuffquant", "Cuffnorm-generated cxb files for all samples.");
#    &fieldlist("unmapped", "Hisat-generated bam files for all samples. Bam files contain only reads that could not be mapped.");
    # end unique to this pipeline
    &space;
    &h1("Acknowledgements");
    &line("Development of this automated analysis pipeline was made possible by support from the University of Minnesota Informatics Institute, in collaboration with the University of Minnesota Genomics Center and the RIS group at the University of Minnesota Supercomputing Institute.");
    &image("umnlogo.gif");
    &space;

#    close $reportFH;
}

####################### helper subs #########################

