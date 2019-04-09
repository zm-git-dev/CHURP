#!/usr/bin/env python
"""Define a class for Sample sheet handling/creation."""

import glob
import sys
import os
import pprint
import re

import GopherPipelines
from GopherPipelines.FileOps import default_files
from GopherPipelines.FileOps import default_dirs
from GopherPipelines.FileOps import dir_funcs
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
        self.sheet_logger.debug(
            'User opts:\n%s', pprint.pformat(self.useropts))
        self.sheet_logger.debug(
            'Default opts:\n%s', pprint.pformat(self.defaultopts))
        self.sheet_logger.debug('Programs for options: %s', self.programs)
        changed = False
        for prog in self.programs:
            if self.useropts[prog]:
                changed = True
                self.finalopts[prog] = self.useropts[prog]
            else:
                self.finalopts[prog] = self.defaultopts[prog]
        self.sheet_logger.debug(
            'Resolved opts:\n%s', pprint.pformat(self.finalopts))
        if changed:
            self.sheet_logger.warning((
                'Be cautious when specifying alternate option strings! '
                'We do not guarantee that they will work. '
                'Always check the manual for the version of the programs that '
                'you are using. This is not an error message.'
                ))
        return

    def _get_fq_paths(self, d):
        """Read through the contents of a FASTQ directory and try to build a
        list of samples from it."""
        # Write a regular expression that will match the parts of the filename
        # that come after the sample name
        samp_re = re.compile(
            r'(_S[0-9]+)?'
            r'(_L00[1-8])?'
            r'(_R(1|2))?_001\.((fq(\.gz)?$)|(fastq(\.gz)?$))')
        # Get all files that look like not-R2 fastq files, make the matching
        # case-insensitive.
        fq_re = re.compile(
            r'^.+[^_R2]_001\.((fq(\.gz)?$)|(fastq(\.gz)?$))',
            flags=re.I)
        cont = os.listdir(d)
        # From the Illumina BaseSpace online documentation, this is what the
        # standard filenames will look like:
        #   SampleName_SX_L00Y_R1_001.fastq.gz
        # X: Sample nummber in samplesheet
        # Y: Lane number
        # R1/2: Fwd/reverse
        # 001: Always 001.
        # This is similar to the UMGC filenames, which are split across lanes
        # and then concatenated.
        # We will iterate through all files in the directory. If they look like
        # R1 fastq files, we will extract the samplename from them. We will
        # then build the R2 filename and ask if that exists in the directory
        for f in cont:
            if re.match(fq_re, f):
                # Extract the samplename from the fastq name
                sn = re.sub(samp_re, '', f)
                # Tack it onto the samples dictionary
                self.samples[sn] = {}
                self.samples[sn]['R1'] = os.path.join(d, f)
                # Look for the R2. This is really dumb-looking but:
                #   Reverse the R1 filename ([::-1])
                #   Replace 1R with 2R, with at most 1 replacement
                #   Reverse it again
                r2 = f[::-1].replace('1R', '2R', 1)[::-1]
                # Extract the samplename from the hypothetical R2 path. If it
                # is different from the R1 samplename, then we have messed up
                # the part of the filename that we shouldn't have - the R2 does
                # not exist for this sammple, and it is single-end
                r2_sn = re.sub(samp_re, '', r2)
                if r2_sn != sn:
                    self.samples[sn]['R2'] = ''
                elif r2 not in cont or r2 == f:
                    self.samples[sn]['R2'] = ''
                elif r2 in cont and r2 != f:
                    self.samples[sn]['R2'] = os.path.join(d, r2)
                else:
                    self.samples[sn]['R2'] = ''
        self.sheet_logger.debug(
            'Found samples:\n%s',
            pprint.pformat(self.samples))
        return

    def write_sheet(self, od, pn, delim):
        """Write the sheet into the given directory."""
        # Define a samplesheet output name
        ssname = default_files.default_samplesheet(pn)
        ssname = os.path.join(od, ssname)
        if os.path.isfile(ssname):
            self.sheet_logger.warning(
                'Samplesheet %s exists. Overwriting', ssname)
        # Write the data into the sheet
        handle = open(ssname, 'w')
        with open(ssname, 'w') as f:
            for sample in sorted(self.final_sheet):
                towrite = [sample] + [
                    self.final_sheet[sample][c]
                    for c
                    in self.column_order]
                handle.write(delim.join(towrite) + '\n')
        # Add a new line
        handle.write('\n')
        # Add the software version data to the bottom of the samplesheet
        handle.write('# Generated by CHURP version ' +
                     GopherPipelines.__version__ + '\n')
        handle.write('# Generated at ' + GopherPipelines.NOW + '\n')
        # And write the major version as a comment for the PBS script to parse
        handle.write('#' + GopherPipelines.__version__.split('.')[0] + '\n')
        handle.flush()
        handle.close()
        return ssname
