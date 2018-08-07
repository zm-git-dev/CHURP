#!/usr/bin/env python
"""Add arguments for the bulk RNAseq pipeline."""

import argparse
from GopherPipelines.FileOps import default_dirs

# Long help messages as constants
EXTRA_OPT_HELP = """Must be passed as a quoted string with = after the option.

Example:
--trimmomatic-opts="-phred64 -threads 4".

This is unfortunately due to a known bug in argparse that prevents proper
handling of options with dashes in them."""


def add_args(ap):
    """Takes an ArgumentParser object, and adds some arguments to it. These will
    be for the bulk RNAseq pipeline. This function returns None; it is called
    for the side-effect of adding arguments to the parser object."""
    # This group is called "required" but it actually isn't set as required by
    # the parser. We will manually check later for the proper combination of
    # arguments that are set.
    ap_req = ap.add_argument_group(
        title='Required Arguments')
    ap_req.add_argument(
        '--fq-folder',
        '-f',
        metavar='<fastq folder>',
        dest='fq_folder',
        required=False,
        help='Directory that contains the FASTQ files.')
    ap_req.add_argument(
        '--hisat2-index',
        '-x',
        metavar='<HISAT2 index base>',
        dest='hisat2_idx',
        required=False,
        help='Basename of the HISAT2 index.')
    ap_req.add_argument(
        '--gtf',
        '-g',
        metavar='<GTF>',
        dest='gtf',
        required=False,
        help='GTF that accompanies the genome for the HISAT2 index.')
    # Make an argument group for optional arguments
    ap_opt = ap.add_argument_group(
        'Optional Arguments')
    ap_opt.add_argument(
        '--help',
        '-h',
        help='Show this help message and exit.',
        action='help')
    ap_opt.add_argument(
        '--verbosity',
        '-v',
        metavar='<loglevel>',
        dest='verbosity',
        help='How much logging output to show. Choose one of "debug," "info," or "warn."',
        choices=['debug', 'info', 'warn'],
        default='warn')
    ap_opt.add_argument(
        '--no-trim',
        dest='no_trim',
        help='If supplied, do not trim reads with Trimmomatic.',
        action='store_true',
        default=False)
    ap_opt.add_argument(
        '--output-dir',
        '-o',
        metavar='<output directory>',
        dest='output_dir',
        help='Output directory. Defaults to global scratch.',
        default=default_dirs.default_output('bulk_rnaseq'))
    ap_opt.add_argument(
        '--threads',
        '-t',
        metavar='<num. threads>',
        dest='threads',
        help='Number of threads to supply to each pipeline step. Default 1.',
        type=int,
        default=1)
    ap_opt.add_argument(
        '--trimmomatic-opts',
        metavar='<trimmomatic options>',
        help='Trimmomatic options. ' + EXTRA_OPT_HELP,
        dest='trimmomatic',
        type=str,
        default='')
    ap_opt.add_argument(
        '--hisat2-opts',
        metavar='<HISAT2 options>',
        help='HISAT2 options. ' + EXTRA_OPT_HELP,
        dest='hisat2',
        type=str,
        default='')
    return
