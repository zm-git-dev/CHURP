#!/usr/bin/env python
"""Define a general Pipeline object that we will use to build the actual
analysis pipelines. This class will have all methods that are common to the
analysis pipelines that we will build: setting dependencies, checking that
programs exist and checking sample lists."""


from GopherPipelines.ArgHandling import set_verbosity


class Pipeline(object):
    """Define the Pipeline object. The following attributes are set:

        self.fq_dir (pathlib.Path): The path to the FASTQ directory
        self.outdir (pathlib.Path): The output directory
        self.required_mods (list): List of required software modules

    The following methods are also defined:

        check_dirs(): Check that directories exist and can be written to.
        prepare_samplesheet(): Initialize a SampleSheet object and populate it
        setup_workdir(): Prepare the working directory
        prepare_qsub(): Build the qsub command line

    """
    def __init__(self, args):
        """Initialize the pipeline object with user-supplied inputs. The
        general pipeline attributes that get set here are:
            - Output directory
            - Program names for options setting
            - User options dictionary
            - Default optiond dictionary
            - Final options dictionary"""
        self.logger = set_verbosity.verb(args['verbosity'], __name__)
        self.logger.debug('Passed args: %s', args)
        self.outdir = args['outdir']
        # These are empty, and will get populated by the sub-class.
        self.programs = []
        self.useropts = {}
        self.defaultopts = {}
        self.finalopts = {}
        return

    def check_dirs(self):
        """Check that the directories exist and are readable and writeable. This
        will raise an error if we cannot find the fastq directory, or the
        output directory cannot be written to."""
        self.logger.info('Checking directories.')
        self.logger.debug('Output dir: ' + str(self.outdir))
        pass

    def resolve_options(self):
        """Read the options dictionary that is set by the user and by the
        sub-class initialization. Basically, any supplied user options will
        clobber the default options. We WILL NOT check that the supplied user
        options are valid. We will try to provide ample warning that users are
        responsible for giving options that will work."""
        self.logger.info('Resolving user and default options.')
        self.logger.debug('User opts: %s', self.useropts)
        self.logger.debug('Default opts: %s', self.defaultopts)
        self.logger.debug('Programs for options: %s', self.programs)
        for prog in self.programs:
            if self.useropts[prog]:
                self.finalopts[prog] = self.useropts[prog]
            else:
                self.finalopts[prog] = self.defaultopts[prog]
        self.logger.debug('Resolved opts: %s', self.finalopts)
        if self.finalopts != self.defaultopts:
            self.logger.warning((
                'Be cautious when specifying alternate option strings! '
                'We do not guarantee that they will work. '
                'Always check the manual for the latest version of the '
                'programs that you are using. This is not an error message.'
                ))
        return

    def prepare_samplesheet(self, ss=None):
        """Read the provided UMGC report sheet and build a list of Sample
        objects that will go downstream to the actual job submission routine.
        Optionally, if 'ss' not provided, we will attempt to build one by
        searching the fastq folder and matching filenames. This is not the
        preferred way, though."""
        self.logger.info('Preparing samplesheet.')
        if not ss:
            self.logger.debug('Attempting to build one from ' + str(self.fq_dir))
        else:
            self.logger.info('Parsing ' + str(ss))
        pass

    def setup_workdir(self):
        """Make the output directory and make sure that the needed materials
        are present."""
        self.logger.info('Setting up working direcotry.')
        pass

    def prepare_qsub(self):
        """Take the pipeline attributes and prepare the qsub call. This drop the
        variables for scratch directory, samplesheet, and index into the qsub
        call."""
        self.logger.info('Setting up qsub command.')
        pass
