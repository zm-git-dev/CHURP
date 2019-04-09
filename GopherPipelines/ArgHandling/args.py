#!/usr/bin/env python
"""Define a function that will parse the arguments. This will eventually be
a huge function, so we will isolate it into this script for ease of development
and maintenance."""

import argparse
import GopherPipelines
from GopherPipelines.ArgHandling import group_template_args
from GopherPipelines.ArgHandling import bulk_rnaseq_args
from GopherPipelines.ArgHandling import sc_rnaseq_args
from GopherPipelines import DieGracefully


# Longer help messages as constants
PIPE_HELP = 'CHURP pipelines for high-throughput sequencing data analysis.'
LIST_HELP = ('List species shorthands for automatically setting alignment '
             'targets and annotation databases.')
GROUP_HELP = ('Make sample metadata CSV templates for statistical tests. These'
              ' templates can be used to make files that designate groups for '
              'differential expression testing, for example.')
BRNASEQ_HELP = ('Bulk RNAseq analysis, including QC, mapping, and expression. '
                'Alignment is performed with HISAT2. Read counts are '
                'generated with featureCounts from the Subread package. '
                'Expression anaylsis is done with edgeR in R.')
SCRNASEQ_HELP = 'Single-cell RNAseq analysis.'


def usage():
    """Print a usage message for the pipeline. This is invoked when there are no
    arguments supplied to the script."""
    msg = """Usage: churp.py <subcommand> <options>

where <subcommand> is the name of the pipeline that is to be run. The specified
<options> will be applied to the operations in the pipeline. Each pipeline has
its own set of options that must be specified. To see a full listing of each
available option for a given pipeline, pass the '--help' option. Alternately,
online help is maintained at the GitHub repository.

Currently, the following subcommands are supported:
    - group_template
    - bulk_rnaseq

For issues, contact help@msi.umn.edu.
Version: {version}
{date}"""
    print(
        msg.format(
            version=GopherPipelines.__version__,
            date=GopherPipelines.__date__))
    return


def check_for_bad(a):
    """Check for bad characters in the args and throw an error if it detects
    any of them. These are characters that let users terminate the current
    command and start another, e.g., a file called 'sample_01; rm -rf ~'"""
    v = vars(a)
    bad_chars = [';', '#', '|', '$(', '<', '>', '`']
    for option in v:
        if not v[option]:
            continue
        else:
            for bc in bad_chars:
                if bc in str(v[option]):
                    DieGracefully.die_gracefully(
                        DieGracefully.NEFARIOUS_CHAR, str(v[option]))
    return


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

    # Species list parser
    list_parser = pipe_parser.add_parser(
        'list_species',
        help=LIST_HELP,
        add_help=False)
    # Group template parser
    group_parser = pipe_parser.add_parser(
        'group_template',
        help=GROUP_HELP,
        add_help=False)
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

    group_template_args.add_args(group_parser)
    bulk_rnaseq_args.add_args(bulk_rnaseq_parser)
    sc_rnaseq_args.add_args(sc_rnaseq_parser)
    pargs = parser.parse_args()
    check_for_bad(pargs)
    return vars(pargs)
