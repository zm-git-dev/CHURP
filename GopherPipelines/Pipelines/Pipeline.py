#!/usr/bin/env python
"""Define a general Pipeline object that we will use to build the actual
analysis pipelines. This class will have all methods that are common to the
analysis pipelines that we will build: setting dependencies, checking that
programs exist and checking sample lists."""


class Pipeline(object):
    """Define the Pipeline object. The following attributes are set:

    The following methods are also defined:

    """
    def __init__(self, fq, out):
        """Initialize the pipeline object with user-supplied inputs. The
        general pipeline attributes that get set here are:
            - path to FASTQ folder
            - Output directory"""
        self.fq_dir = fq
        self.outdir = out
        # We will always want to run FastQC on the data.
        self.required_mods = ['fastqc']
        return
