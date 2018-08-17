#!/usr/bin/env python
"""Define a class for Sample sheet handling/creation."""

import glob
import sys
import os
import datetime
import getpass
import pprint

from GopherPipelines.FileOps import default_dirs
from GopherPipelines.ArgHandling import bulk_rnaseq_args
from GopherPipelines.ArgHandling import sc_rnaseq_args

class Samplesheet(object):
    """Holds the data about samples for the pipelines."""
    # These are going to be defined for every samplesheet:
    #   Programs to merge options for
    #   User supplied options
    #   Default options
    #   Final (merged) options
    #   Samples to process
    #   Order of the columns
    #   Final sheet dictionary
    programs = []
    useropts = {}
    defaultopts = {}
    finalopts = {}
    samples = {}
    column_order = []
    final_sheet = {}

    def _resolve_options(self):
        """Read the options dictionary that is set by the user and by the
        sub-class initialization. Basically, any supplied user options will
        clobber the default options. We WILL NOT check that the supplied user
        options are valid. We will try to provide ample warning that users are
        responsible for giving options that will work."""
        self.sheet_logger.info('Resolving user and default options.')
        self.sheet_logger.debug('User opts:\n%s', pprint.pformat(self.useropts))
        self.sheet_logger.debug('Default opts:\n%s', pprint.pformat(self.defaultopts))
        self.sheet_logger.debug('Programs for options: %s', self.programs)
        changed = False
        for prog in self.programs:
            if self.useropts[prog]:
                changed = True
                self.finalopts[prog] = self.useropts[prog]
            else:
                self.finalopts[prog] = self.defaultopts[prog]
        self.sheet_logger.debug('Resolved opts:\n%s', pprint.pformat(self.finalopts))
        if changed:
            self.sheet_logger.warning((
                'Be cautious when specifying alternate option strings! '
                'We do not guarantee that they will work. '
                'Always check the manual for the version of the programs that '
                'you are using. This is not an error message.'
                ))
        return

    def write_sheet(self, od, pn, delim):
        """Write the sheet into the given directory."""
        # Define a samplesheet output name
        today = datetime.date.isoformat(datetime.datetime.today())
        uname = getpass.getuser()
        ssname = today + '.' + uname + '.' + pn + '.samplesheet.txt'
        ssname = os.path.join(od, ssname)
        if os.path.isfile(ssname):
            newname = ssname.replace('.txt', '.rerun.txt')
            self.sheet_logger.warning('Samplesheet %s exists, making a new copy at %s', ssname, newname)
            ssname = newname
        # Write the data into the sheet
        handle = open(ssname, 'w')
        handle.write(delim.join(['SampleNM'] + self.column_order) + '\n')
        with open(ssname, 'w') as f:
            for sample in sorted(self.final_sheet):
                towrite = [sample] + [self.final_sheet[sample][c] for c in self.column_order]
                handle.write(delim.join(towrite) + '\n')
        handle.flush()
        handle.close()
        return ssname

