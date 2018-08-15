#!/usr/bin/env python
"""Define a sub-class of the Pipeline class for bulk RNAseq analysis."""

import pprint
import os
import glob
import subprocess

from GopherPipelines import DieGracefully
from GopherPipelines.Pipelines import Pipeline
from GopherPipelines.ArgHandling import set_verbosity
from GopherPipelines.FileOps import species_list


class BulkRNAseqPipeline(Pipeline.Pipeline):
    """Sub-class of Pipeline for bulk RNAseq analysis."""

    # Define the pipeline name here
    pipe_name = 'bulk_rnaseq'

    def setup(self, args):
        """Initialize the pipeline object. We will call the general
        Pipeline.__init__() here, as well as set some specific pipeline
        attributes."""
        # Set up the verbosity and logging
        self.pipe_logger = set_verbosity.verb(args['verbosity'], __name__)
        # First validate the arguments. We do this in the subclass because the
        # dependencies of the arguments are pipeline-specific.
        valid_args = self._validate_args(args)
        self.pipe_logger.debug('New BulkRNAseqPipeline instance.')
        self.pipe_logger.debug('Validated args:\n%s', pprint.pformat(valid_args))

        # Set fastq directory here
        self.adapters = valid_args['adapters']
        # This pipeline takes options for trimmomatic and hisat2
        self.programs.extend(['trimmomatic', 'hisat2'])
        # Set the default trimmomatic options here. This is from YZ's scripts
        self.defaultopts['trimmomatic'] = 'ILLUMINACLIP:' + self.adapters + ':4:15:7:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18'
        self.defaultopts['hisat2'] = ''
        # Set the user options here
        self.useropts['trimmomatic'] = valid_args['trimmomatic']
        self.useropts['hisat2'] = valid_args['hisat2']
        # Set the overwrite flag
        if valid_args['overwrite']:
            self.ow = '1'
        else:
            self.ow = '0'

        # Set the paths to the single sample analysis script. This will be
        # submitted to the scheduler. This is a little ugly, but because the
        # package contains the shell script data, it should be OK to define a
        # relative path. __file__ is the currently running script,
        # os.path.realpath() gives the full path, with symlinks resolved. We
        # then rsplit() the string returned by realpath() to get the base dir
        # of the gopher-pipelines scipt. Woof.
        self.single_sample_script = os.path.join(
            os.path.realpath(__file__).rsplit(os.path.sep, 3)[0],
            'PBS',
            'bulk_rnaseq_single_sample.pbs')
        self.summary_script = os.path.join(
            os.path.realpath(__file__).rsplit(os.path.sep, 3)[0],
            'PBS',
            'run_summary_stats.pbs')
        return

    def _validate_args(self, a):
        """Validate arguments for the BulkRNAseqPipeline object. We define it in
        this file because it only really needs to be accessible to this subclass.
        Argument dependencies are pipeline-specific. For the bulk RNAseq analysis
        pipeline:
            - If organism was specified, then the HISAT2 index and GTF do not
              need to be checked
              - Check that the HISAT2 index is complete
        """
        # Check the completeness of the argument dictionary. Either -f or
        # -u must be specified. (-x and -g) or -r must be specified.
        if not a['fq_folder'] and not a['umgc']:
            DieGracefully.die_gracefully(DieGracefully.BRNASEQ_INC_ARGS)
        elif not ((a['hisat2_idx'] and a['gtf']) or a['organism']):
            DieGracefully.die_gracefully(DieGracefully.BRNASEQ_INC_ARGS)
        elif (a['hisat2_idx'] and a['gtf']) and a['organism']:
            DieGracefully.die_gracefully(DieGracefully.BRNASEQ_CONFLICT)
        # Convert all of the paths into absolute paths
        a['outdir'] = os.path.realpath(os.path.expanduser(str(a['outdir'])))
        a['workdir'] = os.path.realpath(os.path.expanduser(str(a['workdir'])))
        a['hisat2_idx'] = os.path.realpath(os.path.expanduser(str(a['hisat2_idx'])))
        # Validate the hisat2 index
        self._validate_hisat_idx(a['hisat2_idx'])
        return a

    def _validate_hisat_idx(self, i):
        """Raise an error if the provided HISAT2 index is not complete -
        all of the [1-8].ht2l? files should be present."""
        # Build glob patterns for the normal and long indices
        norm_idx = i + '.[1-8].ht2'
        long_idx = i + '.[1-8].ht2l'
        # Do the search
        self.pipe_logger.debug('Searching for %s', norm_idx)
        norm_idx_files = glob.glob(norm_idx)
        self.pipe_logger.debug('Found %i index files', len(norm_idx_files))
        # There should be 8 total
        if len(norm_idx_files) == 8:
            return
        else:
            self.pipe_logger.debug('Normal index not found. Searching for long index.')
            long_idx_files = glob.glob(long_idx)
            self.pipe_logger.debug('Found %i long index files', len(long_idx_files))
            if len(long_idx_files) == 8:
                return
            else:
                self.pipe_logger.error('Cound not find HISAT2 index files!')
                DieGracefully.die_gracefully(DieGracefully.BAD_HISAT)
        return

    
    def qsub(self):
        """Write the qsub command. We will need the path to the samplesheet,
        the number of samples in the samplesheet, and the scheduler options that
        were passed as arguments. This is defined in the subclass rather than
        in the main Pipeline class because the exact form of the qsub command
        depends on which pipeline we are running."""
        self.pipe_logger.info('Setting up qsub command.')
        # Set some dummy values here just for testing purposes
        nsamp = 5
        ss_path = '/panfs/roc/scratch/test/samplesheet.txt'
        aln_job_id = '12345'
        self.pipe_logger.debug('Number of samples: %i', nsamp)
        self.pipe_logger.debug('Samplesheet: %s', ss_path)
        # This is the command for aligning and cleaning
        qsub_resources = '"mem=' + str(self.mem) + 'mb'
        qsub_resources += ',nodes=1:ppn=' + str(self.ppn)
        qsub_resources += ',walltime=' + str(self.walltime * 3600) + '"'
        aln_cmd = [
            'qsub',
            '-q',
            'mesabi',
            '-l',
            qsub_resources,
            '-t',
            '1-' + str(nsamp),
            '-v',
            '"SampleSheet='+ss_path+',overwrite=' + self.ow + '"',
            self.single_sample_script]
        # This is the command for counting and normalizing reads
        summary_cmd = [
            'qsub',
            '-q',
            'mesabi',
            '-l',
            qsub_resources,
            '-W',
            'depend=afterok:' + aln_job_id,
            '-v',
            '"base_in='+ss_path+'"',
            self.summary_script]
        self.pipe_logger.debug('qsub:\n%s', ' '.join(aln_cmd))
        self.pipe_logger.debug('qsub:\n%s', ' '.join(summary_cmd))
        return
