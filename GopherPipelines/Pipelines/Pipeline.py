#!/usr/bin/env python
"""Define a general Pipeline object that we will use to build the actual
analysis pipelines. This class will have all methods that are common to the
analysis pipelines that we will build: setting dependencies, checking that
programs exist and checking sample lists."""

import pprint

from GopherPipelines.ArgHandling import set_verbosity
from GopherPipelines.SampleSheet import SampleSheet
from GopherPipelines.FileOps import dir_funcs

# Define exit codes as constants for diagnostic problems
BAD_OUTDIR = 10
BAD_WORKDIR = 11
BAD_RESOURCES = 12


class Pipeline(object):
    """Define the Pipeline object."""

    def __init__(self, args):
        """Initialize the pipeline object with user-supplied inputs. The
        general pipeline attributes that get set here are:
            - Output directory
            - Program names for options setting
            - User options dictionary
            - Default optiond dictionary
            - Final options dictionary
            - Path to single sample PBS script"""
        self.logger = set_verbosity.verb(args['verbosity'], __name__)
        self.logger.debug('Passed args:\n%s', pprint.pformat(args))
        self.outdir = args['outdir']
        self.workdir = args['workdir']
        # These are empty, and will get populated by the sub-class.
        self.programs = []
        self.useropts = {}
        self.defaultopts = {}
        self.finalopts = {}
        self.single_sample_script = ''
        # Initialize the samplesheet here
        self.sheet = SampleSheet.Samplesheet(args)
        # Set the scheduler resources and PBS options here
        self.ppn = args['ppn']
        self.mem = args['mem']
        self.walltime = args['walltime']
        self.overwrite = args['overwrite']
        return

    def _run_checks(self):
        """Run the checks on the directories and options that were passed.
        This happens after initialization of sub-classes from the pipeline."""
        self._check_scheduler()
        self._check_dirs()
        self._resolve_options()
        return

    def _check_scheduler(self):
        """Check that the scheduler resource requests make sense. ppn should be
        between 1 and 24; mem should be between 2000 and 62000; walltime should
        be between 2h and 96h."""
        try:
            assert self.ppn >= 1 and self.ppn <= 24
        except AssertionError as e:
            self.logger.error(
                'PPN value of %i is invalid! Please specify between 1 and 24.',
                self.ppn)
            exit(BAD_RESOURCES)
        try:
            assert self.mem >= 1 and self.mem <= 62000
        except AssertionError as e:
            self.logger.error(
                'Mem value of %i is invalid! Please specify between 1 and 62000.',
                self.mem)
            exit(BAD_RESOURCES)
        try:
            assert self.walltime >= 1 and self.walltime <= 96
        except AssertionError as e:
            self.logger.error(
                'Walltime value of %i is invalid! Please specify between 1 and 96.',
                self.walltime)
            exit(BAD_RESOURCES)
        return

    def _check_dirs(self):
        """Check that the directories exist and are readable and writeable. This
        will raise an error if we cannot find the fastq directory, or the
        output directory cannot be written to."""
        self.logger.info('Checking directories.')
        # We need to check the output directory and the working directory
        self.logger.debug('Checking output directory %s', self.outdir)
        # First,check that it exists
        if dir_funcs.dir_exists(self.outdir, self.logger):
            # Is it empty?
            if dir_funcs.dir_empty(self.outdir, self.logger):
                # And lastly, is it writeable?
                if dir_funcs.dir_writeable(self.outdir, self.logger):
                    # All good!
                    self.logger.debug('Output dir %s is valid', self.outdir)
                    pass
                else:
                    self.logger.error('Output dir %s cannot be written to!', self.outdir)
                    exit(BAD_OUTDIR)
            else:
                self.logger.warning('Output dir %s is not empty!', self.outdir)
        else:
            self.logger.warning('Output dir %s does not exist, making it', self.outdir)
            s = dir_funcs.make_dir(self.outdir, self.logger)
            if not s:
                exit(BAD_OUTDIR)

        # And do the same for the work directory
        self.logger.debug('Checking working directory %s', self.workdir)
        # First,check that it exists
        if dir_funcs.dir_exists(self.workdir, self.logger):
            # Is it empty?
            if dir_funcs.dir_empty(self.workdir, self.logger):
                # And lastly, is it writeable?
                if dir_funcs.dir_writeable(self.workdir, self.logger):
                    # All good!
                    self.logger.debug('Working dir %s is valid', self.workdir)
                    pass
                else:
                    self.logger.error('Working dir %s cannot be written to!', self.workdir)
                    exit(BAD_WORKDIR)
            else:
                self.logger.warning('Working dir %s is not empty!', self.workdir)
        else:
            self.logger.warning('Working dir %s does not exist, making it', self.workdir)
            s = dir_funcs.make_dir(self.workdir, self.logger)
            if not s:
                exit(BAD_WORKDIR)
        return

    def _resolve_options(self):
        """Read the options dictionary that is set by the user and by the
        sub-class initialization. Basically, any supplied user options will
        clobber the default options. We WILL NOT check that the supplied user
        options are valid. We will try to provide ample warning that users are
        responsible for giving options that will work."""
        self.logger.info('Resolving user and default options.')
        self.logger.debug('User opts:\n%s', pprint.pformat(self.useropts))
        self.logger.debug('Default opts:\n%s', pprint.pformat(self.defaultopts))
        self.logger.debug('Programs for options: %s', self.programs)
        for prog in self.programs:
            if self.useropts[prog]:
                self.finalopts[prog] = self.useropts[prog]
            else:
                self.finalopts[prog] = self.defaultopts[prog]
        self.logger.debug('Resolved opts:\n%s', pprint.pformat(self.finalopts))
        if self.finalopts != self.defaultopts:
            self.logger.warning((
                'Be cautious when specifying alternate option strings! '
                'We do not guarantee that they will work. '
                'Always check the manual for the version of the programs that '
                'you are using. This is not an error message.'
                ))
        return

    def prepare_samplesheet(self):
        """Call the samplesheet build method here. The SampleSheet object has
        the fastq directory defined within it, so we do not have to pass any
        other data to it."""
        self.logger.info('Preparing samplesheet.')
        self._run_checks()
        print(self.sheet.fq_dir)
        pass

