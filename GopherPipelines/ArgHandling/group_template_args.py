#!/usr/bin/env python
"""Add arguments for the group template building subcommand."""

import argparse

from GopherPipelines.FileOps import default_dirs
from GopherPipelines.FileOps import default_files


def add_args(ap):
    """Takes an ArgumentParser object, and adds arguments to it. These args
    will be for the group template building subcommand. The function returns
    NoneType; we only call it for the side-effect of adding arguments to the
    parser object."""
    # Try to define a nested subparser. This will get wild
    ap_sub = ap.add_subparsers(
        dest='pipe_group',
        title='Pipelines',
        help='Available pipelines')

    # This is a parser for the bulk RNAseq options.
    brnaseq_def_csv = default_dirs.default_output('bulk_rnaseq') / \
        default_files.default_group_csv('bulk_rnaseq')
    brnaseq_group = ap_sub.add_parser(
        'bulk_rnaseq',
        help='Generate template groups file for bulk RNAseq analysis.',
        add_help=False)
    brnaseq_group_req = brnaseq_group.add_argument_group(
        title='Required arguments')
    brnaseq_group_req.add_argument(
        '--fq-folder',
        '-f',
        metavar='<fastq folder>',
        dest='fq_folder',
        help='Directory that contains the FASTQ files.',
        required=True)
    brnaseq_group_opt = brnaseq_group.add_argument_group(
        title='Optional arguments')
    brnaseq_group_opt.add_argument(
        '--help',
        '-h',
        help='Show this help message and exit.',
        action='help')
    brnaseq_group_opt.add_argument(
        '--verbosity',
        '-v',
        metavar='<loglevel>',
        dest='verbosity',
        help=('How much logging output to show. '
              'Choose one of "debug," "info," or "warn."'),
        choices=['debug', 'info', 'warn'],
        default='warn')
    brnaseq_group_opt.add_argument(
        '--output',
        '-o',
        metavar='<output file>',
        dest='outfile',
        help='Write the CSV to this file. Defaults to ' + str(brnaseq_def_csv),
        default=brnaseq_def_csv)
    brnaseq_group_opt.add_argument(
        '--extra-column',
        '-e',
        metavar='<extra column>',
        action='append',
        default=[],
        help=('Extra column to add to the CSV file. This option may be passed '
              'more than once to add multiple columns. Please note that if '
              'you use this option, then you will need to run a customised '
              'edgeR or DESeq2 analysis to properly handle the experimental '
              'conditions in your metadata CSV. The "SampleName" and "Group" '
              'columns are present by default.'))
    return
