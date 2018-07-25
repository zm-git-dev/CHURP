#!/usr/bin/env python3
"""Gopher-pipelines control script, for running standard analyses of various
types of high-throughput sequencing data analyses. This script is written to
work on the University of Minnesota Supercomputing Institute clusters with
data from the University of Minnesota Genomics Centre. The following analysis
pipelines are supported:
    - Bulk RNAseq
Questions should be directed to help@msi.umn.edu.
Version: 0.0
2017-07-24
"""

# Import standard library modules here
import sys
import os
import argparse
# Import the package here
try:
    from GopherPipelines.ArgHandling import args
except ImportError:
    sys.stderr.write('Error - GopherPipelines package cannot be not found!\n')
    exit(1)

# And check the Python version
try:
    assert sys.version_info.major == 3
except AssertionError:
    sys.stderr.write('Error - this pipeline requires Python 3!\n')
    exit(2)


def brnaseq(args):
    """This function loads the bulk RNAseq pipeline module, and runs through the
    steps for bulk RNAseq analysis."""
    from GopherPipelines.Pipelines import BulkRNAseq
    p = BulkRNAseq.BulkRNAseqPipeline(
        args['fq_folder'],
        args['output_dir'],
        args['no_trim'])
    print(p)
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
        # Choose which pipeline to run here
        if pipe_args['pipeline'] == 'bulk_rnaseq':
            brnaseq(pipe_args)
    return


main()
