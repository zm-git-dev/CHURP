#!/usr/bin/env python
"""Define a general Pipeline object that we will use to build the actual
analysis pipelines. This class will have all methods that are common to the
analysis pipelines that we will build: setting dependencies, checking that
programs exist and checking sample lists."""


class Pipeline(object):
    """Define the Pipeline object. The following attributes are set:

        self.fq_dir (pathlib.Path): The path to the FASTQ directory
        self.outdir (pathlib.Path): The output directory
        self.required_mods (list): List of required software modules

    The following methods are also defined:

        check_dirs(): Check that directories exist and can be written to

    """
    def __init__(self, args):
        """Initialize the pipeline object with user-supplied inputs. The
        general pipeline attributes that get set here are:
            - path to FASTQ folder
            - Output directory"""
        self.fq_dir = args['fq_folder']
        self.outdir = args['output_dir']
        # We will always want to run FastQC on the data.
        self.required_mods = ['fastqc']
        return

    def check_dirs(self):
        """Check that the directories exist and are readable and writeable. This
        will raise an error if we cannot find the fastq directory, or the
        output directory cannot be written to."""
        pass

    def preapre_samplesheet(self, ss=None):
        """Read the provided UMGC report sheet and build a list of Sample
        objects that will go downstream to the actual job submission routine.
        Optionally, if 'ss' not provided, we will attempt to build one by
        searching the fastq folder and matching filenames. This is not the
        preferred way, though."""
        pass

    def setup_workdir(self):
        """Make the output directory and make sure that the needed materials
        are present."""
        pass

    def prepare_qsub(self):
        """Take the pipeline attributes and prepare the qsub call. This drop the
        variables for scratch directory, samplesheet, and index into the qsub
        call."""
        pass
