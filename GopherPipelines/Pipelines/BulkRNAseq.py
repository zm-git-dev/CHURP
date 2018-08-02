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
            'Modules: ' + ','.join(self.required_mods),
            'Methods and Attributes: ' + ','.join(dir(self))])
        return s

    def __init__(self, args):
        """Initialize the pipeline object. We will call the general
        Pipeline.__init__() here, as well as set some specific pipeline
        attributes."""
        # We call the parent class init() here
        super().__init__(args)
        # Add hisat2, samtools, and R to the list of required modules. Also
        # add trimmomatic, conditionally
        self.required_mods.extend(
            ('hisat2', 'samtools', 'R')
            )
        if not args['no_trim']:
            self.required_mods.append('trimmomatic')
        return

    def summarize_fastqc(self):
        """Summarize the FASTQC reports for the run."""
        pass

    def summarize_counts(self):
        """Genearate counts from a list of BAM files and write them into a
        scratch directory. Return the filename to the counts."""
        pass
