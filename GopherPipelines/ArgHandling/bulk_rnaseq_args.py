#!/usr/bin/env python
"""Add arguments for the bulk RNAseq pipeline."""

import argparse
from GopherPipelines.FileOps import default_dirs

# Long help messages as constants
EXTRA_TRIM_HELP = """Must be passed as a quoted string with = after the option.

Example:
--trimmomatic-opts="-phred64 -threads 4".

By default we use the following options:
ILLUMINACLIP:4:15:7:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18.
See the Trimmomatic manual for explanation of these parameters.
"""

EXTRA_HISAT_HELP = """Must be passed as a quoted string with = after the option.

Example:
--hisat2-opts="--phred64 --no-unal"

By default we set the following options: -p (number of threads),
--no-mixed (do not split pairs), and --new-summary (machine-friendly summary).
See the HISAT2 manual for all options available.
"""


def add_args(ap):
    """Takes an ArgumentParser object, and adds some arguments to it. These
    will be for the bulk RNAseq pipeline. This function returns None; it is
    called for the side-effect of adding arguments to the parser object."""
    ap_req = ap.add_argument_group(
        title='Required Arguments')
    ap_req.add_argument(
        '--fq-folder',
        '-f',
        metavar='<fastq folder>',
        dest='fq_folder',
        help='Directory that contains the FASTQ files.')
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
        help=('Organism shorthand for automatically setting a HISAT2 index and'
              ' GTF annotation file. Incompatible with -x and -g.'))
    # Make an argument group for optional arguments
    ap_opt = ap.add_argument_group(
        'Optional Arguments')
    ap_opt.add_argument(
        '--help',
        '-h',
        help='Show this help message and exit.',
        action='help')
    ap_opt.add_argument(
        '--expr-groups',
        '-e',
        help=('CSV from the "group_template" subcommand that lists '
              'experimental groups for each sample. This file must be '
              'hand-edited to include experimental groups to enable '
              'differential expression testing. See the manual for details.'),
        default=None)
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
        '--purge',
        dest='purge',
        help='If supplied, overwrite files from pervious runs.',
        action='store_true',
        default=False)
    ap_opt.add_argument(
        '--min-gene-length',
        '-l',
        dest='mingene',
        help=('Minimum gene length to retain for edgeR analysis, in bp. '
              'Default: 200'),
        type=int,
        default=200)
    ap_opt.add_argument(
        '--unstranded',
        dest='unstranded',
        help=('If supplied, count reads as "unstranded" in featureCounts. '
              'The default RNAseq library prep protocol used by UMGC is '
              'stranded. Default: false'),
        action='store_true',
        default=False)
    ap_opt.add_argument(
        '--adapters',
        '-a',
        metavar='<adapter file.fa>',
        dest='adapters',
        help=('Adapters to use for trimming with Trimmomatic. Defaults to '
              '$TRIMMOMATIC/adapters/all_illumina_adapters.fa. --no-trim '
              'causes these to not be used.'),
        default='$TRIMMOMATIC/adapters/all_illumina_adapters.fa')
    ap_opt.add_argument(
        '--rmdup',
        dest='rmdup',
        help='If supplied, remove duplicates. Default: No duplicate removal.',
        action='store_true',
        default=False)
    ap_opt.add_argument(
        '--no-trim',
        dest='no_trim',
        help='If supplied, do not trim reads. Default: Trim reads.',
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
        help='Trimmomatic options. ' + EXTRA_TRIM_HELP,
        dest='trimmomatic',
        type=str,
        default='')
    ap_opt.add_argument(
        '--hisat2-opts',
        metavar='<HISAT2 options>',
        help='HISAT2 options. ' + EXTRA_HISAT_HELP,
        dest='hisat2',
        type=str,
        default='')

    # Make an argument for scheduler options
    ap_sched = ap.add_argument_group(
        'Scheduler Options')
    ap_sched.add_argument(
        '--group',
        '-A',
        metavar='<MSI group>',
        dest='pbs_group',
        help='MSI group to charge for SU usage. Defaults to primary group.',
        type=str,
        default=None)
    ap_sched.add_argument(
        '--ppn',
        '-p',
        metavar='<procs per node>',
        dest='ppn',
        help='Processors to allocate for each task. Defaults to 6.',
        type=int,
        default=6)
    ap_sched.add_argument(
        '--mem',
        '-m',
        metavar='<mem per job (Mb)>',
        dest='mem',
        help=('Memory, in megabytes, to allocate to each task. '
              'Defaults to 12000 (12GB)'),
        type=int,
        default=12000)
    ap_sched.add_argument(
        '--walltime',
        '-w',
        metavar='<wall clock time (hours)>',
        dest='walltime',
        help=('Walltime, in hours, to allocate to each task. '
              'Defaults to 12 hours.'),
        type=int,
        default=12)
    return
