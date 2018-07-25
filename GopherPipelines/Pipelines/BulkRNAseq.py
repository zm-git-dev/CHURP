#!/usr/bin/env python
"""Define a sub-class of the Pipeline class for bulk RNAseq analysis."""

from GopherPipelines.Pipelines import Pipeline


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
            'Modules: ' + ','.join(self.required_mods)])
        return s

    def __init__(self, fq, out, notrim):
        """Initialize the pipeline object. We will call the general
        Pipeline.__init__() here, as well as set some specific pipeline
        attributes."""
        # We call the parent class init() here
        super().__init__(fq, out)
        # Add hisat2, samtools, and R to the list of required modules. Also
        # add trimmomatic, conditionally
        self.required_mods.append('hisat2')
        self.required_mods.append('samtools')
        self.required_mods.append('R')
        if not notrim:
            self.required_mods.append('trimmomatic')
        return
