#!/usr/bin/env python
"""Define a sub-class of the Pipeline class for bulk RNAseq analysis."""

from GopherPipelines.Pipelines import Pipeline
from GopherPipelines.ArgHandling import set_verbosity


class BulkRNAseqPipeline(Pipeline.Pipeline):
    """Sub-class of Pipeline for bulk RNAseq analysis."""

    # Define the pipeline name here
    pipe_name = 'bulk_rnaseq'

    def __repr__(self):
        """Return a string that can be printed nicely for debugging purposes.
        The user shouldn't need to see this, but developers will probably
        use it."""
        s = '\n'.join([
            'Pipeline: ' + self.pipe_name,
            'FASTQ Dir: ' + str(self.fq_dir),
            'Output Dir: ' + str(self.outdir),
            'User options: ' + str(self.useropts)])
        return s

    def __init__(self, args):
        """Initialize the pipeline object. We will call the general
        Pipeline.__init__() here, as well as set some specific pipeline
        attributes."""
        # We call the parent class init() here
        super().__init__(args)
        # Add hisat2, samtools, and R to the list of required modules. Also
        # add trimmomatic, conditionally
        self.pipe_logger = set_verbosity.verb(args['verbosity'], __name__)
        self.pipe_logger.debug('New BulkRNAseqPipeline instance.')
        # This pipeline takes options for trimmomatic and hisat2
        self.programs.extend(['trimmomatic', 'hisat2'])
        # Set the default trimmomatic options here. This is from YZ's scripts
        self.defaultopts['trimmomatic'] = 'ILLUMINACLIP:adapters.fa:4:15:7:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18'
        self.defaultopts['hisat2'] = ''
        # Set the user options here
        self.useropts['trimmomatic'] = args['trimmomatic']
        self.useropts['hisat2'] = args['hisat2']
        return
