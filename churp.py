#!/usr/bin/env python3
"""CHURP control script, for running standard analyses of various types of
high-throughput sequencing data analyses. This script is written to work on the
University of Minnesota Supercomputing Institute clusters with data from the
University of Minnesota Genomics Centre. The following analysis pipelines are
supported:
    - group_template
    - bulk_rnaseq
Questions should be directed to help@msi.umn.edu.
Version: 0.1.0
2019-03-13
"""

# Check the Python version
try:
    import sys
    assert sys.version_info[0] == 3
except AssertionError:
    sys.stderr.write('Error - this pipeline requires Python 3!\n')
    exit(2)


# Import standard library modules here
import sys
import os
import argparse
# Import the package here
try:
    from GopherPipelines.ArgHandling import args
    from GopherPipelines import DieGracefully
except ImportError:
    sys.stderr.write('Error - GopherPipelines package cannot be not found!\n')
    exit(1)


def sp_dbs(args):
    """This function will read and display a list of available species
    databases that can be used as targets for various pipelines."""
    pass


def expr_group(args):
    """This function will generate experimental group CSV files for input into
    the various pipelines."""

    def brnaseq_group(a):
        """Sub-function for calling the bulk RNAseq group template."""
        from GopherPipelines.ExperimentGroup import BulkRNAseqGroup
        eg = BulkRNAseqGroup.BulkRNAseqGroup(args)
        eg.setup(args)
        eg.write_sheet()
        DieGracefully.die_gracefully(
            DieGracefully.BRNASEQ_GROUP_SUCCESS,
            eg.dest)
        return

    grp_cmd = {
        'bulk_rnaseq': brnaseq_group
        }
    try:
        grp_cmd[args['pipe_group']](args)
    except KeyError:
        # If we get here, then we should get trapped by the no-args help func
        from GopherPipelines.ExperimentGroup import ExpGroup
        e = ExpGroup.ExpGroup(args)
    return


def brnaseq(args):
    """This function loads the bulk RNAseq pipeline module, and runs through
    the steps for bulk RNAseq analysis."""
    from GopherPipelines.Pipelines import BulkRNAseq
    p = BulkRNAseq.BulkRNAseqPipeline(args)
    p.setup(args)
    pipeline_fname, samplesheet_fname, , key_name, qsub_dat = p.qsub()
    if not qsub_dat:
        DieGracefully.die_gracefully(
            DieGracefully.BRNASEQ_SUCCESS,
            pipeline_fname,
            samplesheet_fname,
            key_name)
    elif qsub_dat[2].returncode != 0:
        DieGracefully.die_gracefully(
            DieGracefully.BRNASEQ_SUBMIT_FAIL,
            qsub_dat)
    else:
        DieGracefully.die_gracefully(
            DieGracefully.BRNASEQ_SUBMIT_OK,
            pipeline_fname,
            samplesheet_fname,
            qsub_dat)
    return


def main():
    """The main function. This function is a very high-level function, and
    it should really only have the logic and structure of the pipeline that is
    invoked. The modules and accessory functions within the package should
    do all the actual work."""
    if not sys.argv[1:]:
        args.usage()
        exit(3)
    else:
        pipe_args = args.pipeline_args()
        cmd = {
            'bulk_rnaseq': brnaseq,
            'group_template': expr_group,
            'list_species': sp_dbs
            }
        cmd[pipe_args['pipeline']](pipe_args)
    return


main()
