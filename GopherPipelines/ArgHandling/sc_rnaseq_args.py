#!/usr/bin/env python
"""Add arguments for the single cell RNAseq pipeline."""

import argparse
from GopherPipelines.FileOps import default_dirs

# Long help messages as constants
EXTRA_OPT_HELP = """Must be passed as a quoted string with = after the option.

Example:

--extra-trimmomatic-opts="-phred64 -threads 4".

This is unfortunately due to a known bug in argparse that prevents proper
handling of options with dashes in them."""


def add_args(ap):
    """Takes an ArgumentParser object, and adds some arguments to it. These will
    be for the single cell RNAseq pipeline."""
    ap_req = ap.add_argument_group(
        'Required Arguments')
    ap_req.add_argument(
        '--fq-folder',
        '-f',
        required=True,
        help='Directory that contains the FASTQ files.')
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
        help='Minimum level of logging to show. Default INFO.',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO')
    ap_opt.add_argument(
        '--trim',
        help='If supplied, trim reads with Trimmomatic.',
        action='store_true',
        default=False)
    ap_opt.add_argument(
        '--output-dir',
        '-o',
        help='Output directory. Defaults to global scratch.',
        default=default_dirs.default_output('sc_rnaseq'))
    ap_opt.add_argument(
        '--threads',
        '-t',
        help='Number of threads to supply to each pipeline step. Default 1',
        type=int,
        default=1)
    ap_opt.add_argument(
        '--extra-trimmomatic-opts',
        help='Extra Trimmomatic options. ' + EXTRA_OPT_HELP,
        type=str,
        default='')
    return
