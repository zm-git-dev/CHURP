#!/usr/bin/env python
"""Define a class for Sample sheet handling/creation."""

import glob, sys, os
from GopherPipelines.FileOps import default_dirs
from GopherPipelines.ArgHandling import bulk_rnaseq_args
from GopherPipelines.ArgHandling import sc_rnaseq_args

class Samplesheet(object):
    """Holds the data about samples for the pipelines."""
    # These are going to be defined for every samplesheet
    useropts = {}
    defaultopts = {}
    finalopts = {}
    samples = {}

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
        for prog in self.programs:
            if self.useropts[prog]:
                self.finalopts[prog] = self.useropts[prog]
            else:
                self.finalopts[prog] = self.defaultopts[prog]
        self.sheet_logger.debug('Resolved opts:\n%s', pprint.pformat(self.finalopts))
        if self.finalopts != self.defaultopts:
            self.sheet_logger.warning((
                'Be cautious when specifying alternate option strings! '
                'We do not guarantee that they will work. '
                'Always check the manual for the version of the programs that '
                'you are using. This is not an error message.'
                ))
        return

    def parse_sheet(ss):
        """Method for parsing a provided sample sheet, i.e., if ss in Pipelines.prepare_samplesheet()"""
        with open(ss) as sheet:
            for line in sheet:
                add_relevant_stuff_to_dict

    def create_sheet(ss):
        """Method for creating a sample sheet from scratch given the UMGC-provided data"""
        with open(UMGC_spreadsheetloc) as sheet:
            add_stuff_to_dict
