#!/usr/bin/env python
"""A class to handle the generation of experimental groups for each sample.
This may seem like overkill for now, but eventually we may want to generate
complicated metadata sheets or do some more involved statistical analyses on
groups of samples."""

import re
import os
import pprint

from GopherPipelines import DieGracefully
from GopherPipelines.ArgHandling import set_verbosity
from GopherPipelines.FileOps import dir_funcs


class ExpGroup(object):
    """The class that will hold the experimetnal group data, as well as write
    the sheet to the output directory."""

    def __init__(self, args):
        """Initialize the object."""
        # Validate the arguments
        self._group_help(args)
        self.group_logger = set_verbosity.verb(args['verbosity'], __name__)
        self.samples = {}
        self.columns = ['SampleName', 'Group']
        self.dest = os.path.realpath(os.path.expanduser(args['outfile']))
        self._prepare_output()
        return

    def _group_help(self, args):
        """Do a simple check to see if the pipe_group argument is set. If it
        isn't, print a nice help message."""
        if not args['pipe_group']:
            DieGracefully.die_gracefully(DieGracefully.GROUP_NO_PIPE)
        return

    def _prepare_output(self):
        """Make sure the parent directory of the output directory exists and
        can be written into. Make directories that we have permission to
        make."""
        # Get the dirname of the output file to make sure that the directory
        # exists and can be written to
        par = os.path.dirname(self.dest)
        self.group_logger.info('Checking directories.')
        # We need to check the output directory and the working directory
        self.group_logger.debug('Checking output directory %s', par)
        # First,check that it exists
        if dir_funcs.dir_exists(par, self.group_logger):
            # And lastly, is it writeable?
            if dir_funcs.dir_writeable(par, self.group_logger):
                # All good!
                self.group_logger.debug(
                    'Output dir %s is valid', par)
                pass
            else:
                self.group_logger.error(
                    'Output dir %s cannot be written to!', par)
                DieGracefully.die_gracefully(DieGracefully.BAD_OUTDIR)
        else:
            self.group_logger.warning(
                'Output dir %s does not exist, making it', par)
            s = dir_funcs.make_dir(par, self.group_logger)
            if not s:
                DieGracefully.die_gracefully(DieGracefully.BAD_OUTDIR)
        pass

    def _validate_fastq_folder(self, d):
        """Raise an error if the FASTQ directory does not exist or does not
        have any FASTQ files."""
        try:
            contents = os.listdir(d)
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.BAD_FASTQ)
        # Check if there is at least one file ending in a standard fastq suffix
        fq_pat = re.compile(r'^.+((.fq(.gz)?$)|(.fastq(.gz)?$))')
        has_fastq = False
        for f in contents:
            if re.match(fq_pat, f):
                has_fastq = True
                break
        if has_fastq:
            return
        else:
            DieGracefully.die_gracefully(DieGracefully.EMPTY_FASTQ)
        return

    def _get_sample_names(self, d):
        """Read the contents of the supplied FASTQ directory and parse out the
        sample names."""
        # regular expression to slice out the samplename from the read name
        samp_re = re.compile(
            r'(_S[0-9]+)?(_L00[1-8])?'
            r'(_R(1|2))?'
            r'_001\.((fq(\.gz)?$)|(fastq(\.gz)?$))',
            flags=re.I)
        # regular expression to get files that look like not-R2 FASTQ files
        fq_re = re.compile(
            r'^.+[^_R2]_001\.((fq(\.gz)?$)|(fastq(\.gz)?$))',
            flags=re.I)
        cont = os.listdir(d)
        sd = {}
        for f in cont:
            if re.match(fq_re, f):
                sn = re.sub(samp_re, '', f)
                sd[sn] = {}
                self.group_logger.debug('Found sample %s', sn)
        return sd

    def _build_groups(self, dummy='NULL'):
        """Populate the sample dictionary with dummy values."""
        for sn in self.samples:
            for c in self.columns:
                if c == 'SampleName':
                    self.samples[sn][c] = sn
                else:
                    self.samples[sn][c] = dummy
        self.group_logger.debug(
            'Sample groups:\n%s', pprint.pformat(self.samples))
        return

    def write_sheet(self):
        """Write the sheet to the output file."""
        if os.path.isfile(self.dest):
            self.group_logger.warning(
                'Groups file %s exists! Overwriting!', self.dest)
        fh = open(self.dest, 'w')
        fh.write(','.join(self.columns) + '\n')
        for sn in sorted(self.samples):
            towrite = ','.join(
                [self.samples[sn][c] for c in self.columns]) + '\n'
            fh.write(towrite)
        fh.flush()
        fh.close()
        return
