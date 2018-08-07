#!/usr/bin/env python
"""Define a function that will parse the arguments. This will eventually be
a huge function, so we will isolate it into this script for ease of development
and maintenance."""

import argparse
from GopherPipelines.ArgHandling import bulk_rnaseq_args
from GopherPipelines.ArgHandling import sc_rnaseq_args


# Longer help messages as constants
PIPE_HELP = 'Gopher-pipelines for high-throughput sequencing data analysis.'
BRNASEQ_HELP = 'Bulk RNAseq analysis, including QC, mapping, and expression.'
SCRNASEQ_HELP = 'Single-cell RNAseq analysis.'


def usage():
    """Print a usage message for the pipeline. This is invoked when there are no
    arguments supplied to the script."""
    msg = """Usage: gopher-pipelines.py <subcommand> <options>

where <subcommand> is the name of the pipeline that is to be run. The specified
<options> will be applied to the operations in the pipeline. Each pipeline has
its own set of options that must be specified. To see a full listing of each
available option for a given pipeline, pass the '--help' option. Alternately,
online help is maintained at the GitHub repository.

Currently, the following subcommands are supported:
    - bulk_rnaseq

For issues, contact help@msi.umn.edu.
Version: 0.0
2018-07-24"""
    print(msg)
    return


def validate(a):
    """Validate the arguments. This function will run checks to make sure that
    valid combinations of arguments have been passed. For now, it will check
    that either a samplesheet, or the group of fq/hisat2/gtf were provided."""
    return vars(a)


def pipeline_args():
    """Parse arguments for the pipeline as a whole. This will dispatch the
    arguments to other functions that parse pipeline-specific arguments. This
    is a bit crazy for now, but it should be pretty flexible for adding new
    pipelines and easier to maintain in the long run."""
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
    # sc_rnaseq parser
    sc_rnaseq_parser = pipe_parser.add_parser(
        'sc_rnaseq',
        help=SCRNASEQ_HELP,
        add_help=False)

    bulk_rnaseq_args.add_args(bulk_rnaseq_parser)
    sc_rnaseq_args.add_args(sc_rnaseq_parser)
    pargs = parser.parse_args()
    v_args = validate(pargs)
    return v_args
