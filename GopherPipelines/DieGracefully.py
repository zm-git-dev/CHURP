#!/usr/bin/env python
"""Defines some exit codes and error messages to write to the terminal when
an error is encountered. This will hopefully be more helpful to users than a
Python exception or the debug console messages."""

import sys
import os

# The error codes are constants
BAD_OUTDIR = 10
BAD_WORKDIR = 11
BAD_RESOURCES = 12
BAD_FASTQ = 13
EMPTY_FASTQ = 14
BRNASEQ_INC_ARGS = 20
BRNASEQ_CONFLICT = 21
BAD_HISAT = 22
BAD_GTF = 23
BAD_ADAPT = 24
BRNASEQ_BAD_GPS = 25
BRNASEQ_NO_SAMP_GPS = 26
BRNASEQ_SUCCESS = 27
GROUP_NO_PIPE = 50
GROUP_BAD_COL = 51
BRNASEQ_GROUP_SUCCESS = 52

# Define a series of functions that just write messages to the terminal. We
# will call these with a main error handling function at the end.


def general_error():
    """This will be a generic error message for when something goes wrong that
    we can't identify."""
    msg = """
gopher-pipelines has caught an unidentified error. Please send this error
message, the command you typed, and debugging output to the MSI help desk
(help@msi.umn.edu).\n"""
    sys.stderr.write(msg)
    return


def bad_outdir():
    """Call this function when the user supplies a bad output directory."""
    msg = """
The output directory you have supplied is not suitable. Either you do not have
permissions to write into it, or the disk is full. If you have verified that you
can write into the directory, please contact the help desk (help@msi.umn.edu)
with this error and the debugging output.\n"""
    sys.stderr.write(msg)
    return


def bad_workdir():
    """Call this function when the user supplies a bad working directory."""
    msg = """
The working directory you have supplied is not suitable. Either you do not have
permissions to write into it, or the disk is full. If you have verified that you
can write into the direcotry, please contact the help desk (help@msi.umn.edu)
with this error and the debugging output.\n"""
    sys.stderr.write(msg)
    return


def bad_resources():
    """Call this function when the user supplies illegal PBS resources."""
    msg = """
The resouces you have specified are out of allowable bounds for our system. The
number of processors per node (PPN) should be an integer between 1 and 24. The
allocated memory should be specified in megabytes (MB) as an integer between 1
and 62000. The walltime should be specified in hours as an integer between 1 and
96.\n"""
    sys.stderr.write(msg)
    return


def brnaseq_success(pipe_script, samplesheet):
    """Call this function when the bulk RNAseq pipeline finishes successfully.
    It will include the pipeline script and the samplesheet paths in the
    output message."""
    msg = """
Samplesheet and pipeline script generation complete! Their paths are given
below:

Pipeline script: {pn}
Samplesheet: {ss}

Verify the information in the samplesheet, and run the pipeline script with
bash while logged into Mesabi. You will recieve email notifications of job
start/completion/error at your UMN X500 email address. If you need to submit
an error repor, please contact help[at]msi.umn.edu. Please include the
samplesheet, pipeline script, and the error message with your report.\n"""
    sys.stderr.write(msg.format(pn=pipe_script, ss=samplesheet))
    return


def brnaseq_inc():
    """Call this function when the bulk_rnaseq arguments dictionary is
    incomplete."""
    msg = """Error!

You did not specify sufficient options to run the bulk_rnaseq subcommand of
gopher-pipelines. You must specify a FASTQ directory (-f). Additionally, you
must either specify a path to a HISAT2 index (-x) and GTF (-g), or an organism
name (-r). If you are building a group template file, you need only specify a
FASTQ directory. Please fix your command line and re-run.\n"""
    sys.stderr.write(msg)
    return


def brnaseq_conflict():
    """Call this function when there are conflicting arguments passed to the
    bulk_rnaseq pipeline."""
    msg = """Error!

You have specified conflicting options to the bulk_rnaseq subcommand. HISAT2
index (-x) and GTF (-g) are both incompatible with organism (-r). Pass the help
option (-h) to see all available options. Please fix your command line and
re-run.\n"""
    sys.stderr.write(msg)
    return


def bad_hisat2():
    """Call this function when the user supplies a bad HISAT2 index."""
    msg = """
The HISAT2 index that you have supplied is not valid. Please give the path to
the base of the HISAT2 index (without the .*.ht2 extension). It should be given
in the same way as would be given to HISTA2 directly. See the HISAT2 manual for
information on building an index from a FASTA file.\n"""
    sys.stderr.write(msg)
    return


def bad_fastq():
    """Call this fucntion when the user supplies a bad FASTQ folder."""
    msg = """
The FASTQ directory that you have supplied either does not exist or cannot be
read. If you have verified that you can read the directory contents, please
contact the help desk (help@msi.umn.edu) with tihs error message and the
debugging output.\n"""
    sys.stderr.write(msg)
    return


def empty_fastq():
    """Call this fucntion when the user supplies an empty FASTQ folder."""
    msg = """
The FASTQ directory that you have supplied does not appear to contain any FASTQ
or gzipped FASTQ files. Please ensure that the files in the directory have names
that end in one of the following: .fastq, .fastq.gz, .fq., .fq.gz
(case sensitive).\n"""
    sys.stderr.write(msg)
    return


def bad_gtf():
    """Call this function when the user supplies a GTF that does not exist or
    cannot be read."""
    msg = """
The GTF you have supplied does not exist, or cannot be read. Please check your
path for any typos and that any special characters are properly quoted or
escaped and try again.\n"""
    sys.stderr.write(msg)
    return


def bad_adapter():
    """Call thsi function when the user supplies an adapters file that does not
    exist or cannot be read."""
    msg = """
The adapters file that you have supplied either does not exist or cannot be
read. Please check your path for any typos and that any special characters are
properly quoted or escaped, and try again.\n"""
    sys.stderr.write(msg)
    return


def brnaseq_bad_groups():
    """Call this function when a user does passes an experimental groups file
    that does not exist, or is not readable."""
    msg = """
The groups template CSV file that you have provided does not exist, or cannot
be read. Please check your path for any typos and that any spaces or special
characters are properly quoted or escaped.\n"""
    sys.stderr.write(msg)
    return


def brnaseq_no_sample_groups():
    """Call this function when a user supplies a groups CSV that does not have
    any information for the samples in the FASTQ directory."""
    msg = """
The groups CSV file that you supplied does not contain any information about
the samples that were found in the FASTQ directory. Check that you supplied the
correct directory with -f, and that your sample names match exactly between the
CSV and the FASTQ directory. Use the template from --make-groups-template to
see the exact sample names that gopher-pipelines is using to match samples to
groups.\n"""
    sys.stderr.write(msg)
    return


def group_no_pipe():
    """Call this function when a user runs the group_template pipeline, but
    does not supply a pipeline for which to build a group template."""
    msg = """Usage: gopher-pipelines.py group_template <pipeline> <options>

This subcommand allows you to build templates for experimental metadata files.
These files specify experimental conditions or treatment groups for comparisons
in downstream analyses. The files from this subcommand, for example, can be
used to perform differential gene expression testing among treatment groups.
Takes a pipeline as an argument, followed by modifiers to control the output.
See the pipeline-specific help (-h) or the user manual for details.

Currently available pipelines:
    - bulk_rnaseq\n"""
    sys.stderr.write(msg)
    return


def group_bad_col():
    """Call this function when a user supplies an invalid column name (contains
    a comma) to the group_template subcommand."""
    msg = """
You have specified an invalid column name. Column names cannot contain commas,
because the experimental group file is a comma-separated values file.\n"""
    sys.stderr.write(msg)
    return


def brnaseq_group_success(fname):
    """Call this function when the group_template command finishes properly.
    It will tell the user where to file has been written, and what they need
    to do next."""
    msg = """
Your experimental groups template has been written to {path}

All sample groups and any additional columns have specified have been filled
with NULL. Please edit the file and write in the correct values for your
dataset. Any samples with a "Group" value of NULL will not be used in
downstream differential expression analysis. When you have edited the file to
your liking, supply its path to the "bulk_rnaseq" pipeline with the -e
option to enable group testing.\n"""
    sys.stderr.write(msg.format(path=fname))
    return


def die_gracefully(e, *args):
    """Print user-friendly error messages and exit."""
    err_dict = {
        BAD_OUTDIR: bad_outdir,
        BAD_WORKDIR: bad_workdir,
        BAD_RESOURCES: bad_resources,
        BAD_FASTQ: bad_fastq,
        EMPTY_FASTQ: empty_fastq,
        BAD_HISAT: bad_hisat2,
        BRNASEQ_INC_ARGS: brnaseq_inc,
        BRNASEQ_CONFLICT: brnaseq_conflict,
        BAD_GTF: bad_gtf,
        BAD_ADAPT: bad_adapter,
        BRNASEQ_BAD_GPS: brnaseq_bad_groups,
        BRNASEQ_NO_SAMP_GPS: brnaseq_no_sample_groups,
        GROUP_NO_PIPE: group_no_pipe,
        GROUP_BAD_COL: group_bad_col,
        BRNASEQ_GROUP_SUCCESS: brnaseq_group_success,
        BRNASEQ_SUCCESS: brnaseq_success
        }
    try:
        err_dict[e](*args)
    except KeyError:
        general_error()
    exit(e)
    return
