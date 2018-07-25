#!/usr/bin/env python
"""Define a function that will parse the arguments. This will eventually be
a huge function, so we will isolate it into this script for ease of development
and maintenance."""

import argparse
from GopherPipelines.FileOps import default_dirs


# Longer help messages as constants
PIPE_HELP = 'Gopher-pipelines for high-throughput sequencing data analysis.'
BRNASEQ_HELP = 'Bulk RNAseq analysis, including QC, mapping, and expression.'
EXTRA_OPT_HELP = """Must be passed as a quoted string with = after the option.

Example:

--extra-trimmomatic-opts="-phred64 -threads 4".

This is unfortunately due to a known bug in argparse that prevents proper
handling of options with dashes in them."""


def usage():
    """Print a usage message for the pipeline. This is invoked when there are no
    arguments supplied to the script."""
    msg = """Usage: gopher-pipelines.py <subcommand> <options>

where <subcommand> is the name of the pipeline that is to be run. The specified
<options> will be applied to the operations in the pipeline. Each pipeline has
its own set of options that must be specified. To see a full listing of each
available option for a given pipeline, pass the '--help' option. Alternately,
online help is maintained at the GitHub repository.

Currently, the following pipelines are supported:
    - bulk_rnaseq

For issues, contact help@msi.umn.edu.
Version: 0.0
2018-07-24"""
    print(msg)
    return


def parse_arguments():
    """use the argparse library to parse the arguments to the script. We will
    implement a series of subcommands, one for each pipeline. Each pipeline
    will have its own set of arguments."""
    parser = argparse.ArgumentParser(
        description=PIPE_HELP,
        add_help=True)
    # Add a sub-parser for the subcommands. These will be the pipelines.
    pipe_parser = parser.add_subparsers(
        dest='pipeline',
        title='Pipelines',
        help='Available pipelines')

    # This is the bulk_rnaseq parser
    bulk_rnaseq_parser = pipe_parser.add_parser(
        'bulk_rnaseq',
        help=BRNASEQ_HELP,
        add_help=False)
    # Next, we want to make argument groups to override the default names of
    # "optional" and "positional"
    bulk_rnaseq_parser_req = bulk_rnaseq_parser.add_argument_group(
        'Required Arguments')
    bulk_rnaseq_parser_req.add_argument(
        '--fq-folder',
        '-f',
        required=True,
        help='Directory that contains the FASTQ files.')
    bulk_rnaseq_parser_req.add_argument(
        '--hisat2-index',
        '-x',
        required=True,
        help='Basename of the HISAT2 index.')
    bulk_rnaseq_parser_req.add_argument(
        '--gtf',
        '-g',
        required=True,
        help='GTF that accompanies the genome for the HISAT2 index.')
    # Make an argument group for optional arguments
    bulk_rnaseq_parser_opt = bulk_rnaseq_parser.add_argument_group(
        'Optional Arguments')
    bulk_rnaseq_parser_opt.add_argument(
        '--help',
        '-h',
        help='Show this help message and exit.',
        action='help')
    bulk_rnaseq_parser_opt.add_argument(
        '--verbosity',
        '-v',
        help='Minimum level of logging to show. Default INFO.',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO')
    bulk_rnaseq_parser_opt.add_argument(
        '--trim',
        help='If supplied, trim reads with Trimmomatic.',
        action='store_true',
        default=False)
    bulk_rnaseq_parser_opt.add_argument(
        '--output-dir',
        '-o',
        help='Output directory. Defaults to global scratch.',
        default=default_dirs.default_output('bulk_rnaseq'))
    bulk_rnaseq_parser_opt.add_argument(
        '--threads',
        '-t',
        help='Number of threads to supply to each pipeline step. Default 1',
        type=int,
        default=1)
    bulk_rnaseq_parser_opt.add_argument(
        '--extra-trimmomatic-opts',
        help='Extra Trimmomatic options. ' + EXTRA_OPT_HELP,
        type=str,
        default='')
    bulk_rnaseq_parser_opt.add_argument(
        '--extra-hisat2-opts',
        help='Extra HISAT2 options. ' + EXTRA_OPT_HELP,
        type=str,
        default='')

    args = parser.parse_args()
    # argparse uses a weird Namespace object - we use vars() to turn it into a
    # more familiar dictionary object.
    return vars(args)
