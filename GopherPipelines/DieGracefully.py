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


def brnaseq_inc():
    """Call this function when the bulk_rnaseq arguments dictionary is
    incomplete."""
    msg = """Error!

You did not specify sufficient options to run the bulk_rnaseq subcommand of
gopher-pipelines. You must specify a FASTQ directory (-f). Additionally, you
must either specify a path to a HISAT2 index (-x) and GTF (-g), or an organism
name (-r). Please fix your command line and re-run.\n"""
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


def die_gracefully(e):
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
        BAD_ADAPT: bad_adapter
        }
    try:
        err_dict[e]()
    except KeyError:
        general_error()
    exit(e)
    return
