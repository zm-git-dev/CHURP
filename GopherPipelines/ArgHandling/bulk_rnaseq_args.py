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
    ap_req = ap.add_argument_group(
        title='Required Arguments')
    ap_req_ss = ap_req.add_mutually_exclusive_group()
    ap_req_ss.add_argument(
        '--umgc-sheet',
        '-u',
        metavar='<UMGC sheet.xlsx>',
        dest='umgc',
        help=('UMGC sheet giving R1 and R2 paths. This must come directly from '
              'the UMGC facility as an Excel file. We do not support parsing '
              'any files that have modifications from the standard format. '
              'Incompatible with -f.')
        )
    ap_req_ss.add_argument(
        '--fq-folder',
        '-f',
        metavar='<fastq folder>',
        dest='fq_folder',
        help=('Directory that contains the FASTQ files. If this is provided, '
              'we will attempt to build the list of R1 and R2 read paths. '
              'Incompatible with -u.'))
    ap_req.add_argument(
        '--hisat2-index',
        '-x',
        metavar='<HISAT2 index base>',
        dest='hisat2_idx',
        required=False,
        help='Basename of the HISAT2 index. Incompatible with -r.')
    ap_req.add_argument(
        '--gtf',
        '-g',
        metavar='<GTF.gtf>',
        dest='gtf',
        required=False,
        help=('GTF that accompanies the genome for the HISAT2 index. '
              'Incompatible with -r.'))
    ap_req.add_argument(
        '--organism',
        '-r',
        metavar='<organism>',
        dest='organism',
        required=False,
        help=('Organism shorthand for automatically setting a HISAT2 index and '
              'GTF annotation file. Pass --list-species to see a full list. '
              'Incompatible with -x and -g.'))
    # Make an argument group for optional arguments
    ap_opt = ap.add_argument_group(
        'Optional Arguments')
    ap_opt.add_argument(
        '--help',
        '-h',
        help='Show this help message and exit.',
        action='help')
    ap_opt.add_argument(
        '--list-species',
        help='List species databases that can be passed to -r.',
        action='store_true',
        default=False)
    ap_opt.add_argument(
        '--verbosity',
        '-v',
        metavar='<loglevel>',
        dest='verbosity',
        help=('How much logging output to show. '
              'Choose one of "debug," "info," or "warn."'),
        choices=['debug', 'info', 'warn'],
        default='warn')
    ap_opt.add_argument(
        '--overwrite',
        dest='overwrite',
        help='Overwrite previous runs?',
        action='store_true',
        default=False)
    ap_opt.add_argument(
        '--adapters',
        '-a',
        metavar='<adapter file.fa>',
        dest='adapters',
        help=('Adapters to use for trimming with Trimmomatic. Defaults to '
              '/home/msistaff/public/all_adapter.fa. --no-trim causes these '
              'to not be used.'),
        default='/home/msistaff/public/all_adapter.fa')
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
        dest='outdir',
        help='Output directory. Defaults to global scratch.',
        default=default_dirs.default_output('bulk_rnaseq'))
    ap_opt.add_argument(
        '--working-dir',
        '-d',
        metavar='<working directory>',
        dest='workdir',
        help='Working directory. Defaults to global scratch.',
        default=default_dirs.default_workdir('bulk_rnaseq'))
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

    # Make an argument for scheduler options
    ap_sched = ap.add_argument_group(
        'Scheduler Options')
    ap_sched.add_argument(
        '--ppn',
        '-p',
        metavar='<procs per node>',
        dest='ppn',
        help='Processors to allocate for each task. Defaults to 8.',
        type=int,
        default=8)
    ap_sched.add_argument(
        '--mem',
        '-m',
        metavar='<mem per job (Mb)>',
        dest='mem',
        help='Memory, in megabytes, to allocate to each task. Defaults to 24000 (24GB)',
        type=int,
        default=24000)
    ap_sched.add_argument(
        '--walltime',
        '-w',
        metavar='<wall clock time (hours)>',
        dest='walltime',
        help='Walltime, in hours, to allocate to each task. Defaults to 12 hours.',
        type=int,
        default=12)
    return
