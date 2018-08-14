#!/usr/bin/env python
"""Define a sub-class of the Pipeline class for bulk RNAseq analysis."""

import pprint
import os
import subprocess

from GopherPipelines.Pipelines import Pipeline
from GopherPipelines.ArgHandling import set_verbosity
from GopherPipelines.FileOps import species_list


class BulkRNAseqPipeline(Pipeline.Pipeline):
    """Sub-class of Pipeline for bulk RNAseq analysis."""

    # Define the pipeline name here
    pipe_name = 'bulk_rnaseq'

    def __init__(self, args):
        """Initialize the pipeline object. We will call the general
        Pipeline.__init__() here, as well as set some specific pipeline
        attributes."""
        # First validate the arguments. We do this in the subclass because the
        # dependencies of the arguments are pipeline-specific.
        valid_args = self.validate_args(args)
        # We call the parent class init() here
        super().__init__(valid_args)
        # Set up the verbosity and logging
        self.pipe_logger = set_verbosity.verb(valid_args['verbosity'], __name__)
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

    def validate_args(self, a):
        """Validate arguments for the BulkRNAseqPipeline object. We define it in
        this file because it only really needs to be accessible to this subclass.
        Argument dependencies are pipeline-specific. For the bulk RNAseq analysis
        pipeline:
            - If organism was specified, then the HISAT2 index and GTF do not
              need to be checked
        """
        return a

    def validate_hisat_idx(self):
        """Raise an error if the provided HISAT2 index is not complete -
        all of the [1-8].ht2l? files should be present."""
        pass

    
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
        qsub_resources = 'mem=' + str(self.mem) + 'mb'
        qsub_resources += ',nodes=1:ppn=' + str(self.ppn)
        qsub_resources += ',walltime=' + str(self.walltime * 3600) + '"'
        aln_cmd = [
            'qsub',
            '-l',
            qsub_resources,
            '-t',
            '1-' + str(nsamp),
            '-v',
            '"SampleSheet='+ss_path+'"',
            self.single_sample_script]
        # This is the command for counting and normalizing reads
        summary_cmd = [
            'qsub',
            '-l',
            qsub_resources,
            '-W',
            'depend=afterok:' + aln_job_id,
            '-v',
            '"base_in='+ss_path+'"',
            self.summary_script]
        self.pipe_logger.debug('qsub: %s', aln_cmd)
        self.pipe_logger.debug('qsub: %s', summary_cmd)
        return
