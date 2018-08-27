#!/usr/bin/env python
"""Define a sub-class of samplesheet that holds the data for bulk RNAseq
samples."""

import os
import re
import pprint

from GopherPipelines.SampleSheet import SampleSheet
from GopherPipelines.ArgHandling import set_verbosity


class BulkRNASeqSampleSheet(SampleSheet.Samplesheet):
    """A sub-class of SampleSheet that holds information for the bulk RNAseq
    pipelines."""

    def __init__(self, args):
        """Initialize the bulk RNAseq samplesheet."""
        # Set up a logger
        self.sheet_logger = set_verbosity.verb(args['verbosity'], __name__)
        # This pipeline takes options for trimmomatic and hisat2
        self.programs.extend(['trimmomatic', 'hisat2'])
        # Set the default trimmomatic options here. This is from YZ's scripts
        self.defaultopts['trimmomatic'] = ''.join([
            'ILLUMINACLIP:',
            args['adapters'],
            ':4:15:7:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18']
            )
        self.defaultopts['hisat2'] = ''
        # Set the user options here
        self.useropts['trimmomatic'] = args['trimmomatic']
        self.useropts['hisat2'] = args['hisat2']
        if args['no_trim']:
            self.useropts['trim'] = 'no'
        else:
            self.useropts['trim'] = 'yes'
        if args['stranded']:
            self.useropts['stranded'] = 'yes'
        else:
            self.useropts['stranded'] = 'no'
        self.useropts['gtf'] = args['gtf']
        self.useropts['hisat2_idx'] = args['hisat2_idx']
        self.useropts['hisat2_threads'] = '-p ' + str(args['ppn']) + ' '
        # Set the column order to be the columns of the sample sheet. This will
        # eventually become the header of the sheet.
        self.column_order.extend([
            'FastqR1files',
            'FastqR2file',
            'OutputDir',
            'WorkingDir',
            'TRIM',
            'trimmomaticOpts',
            'Hisat2index',
            'Hisat2Options',
            'Stranded',
            'AnnotationGTF'])
        self._get_fq_paths(args['fq_folder'])
        self._resolve_options()
        return

    def _get_fq_paths(self, d):
        """Read through the contents of a FASTQ directory and try to build a
        list of samples from it."""
        # Write a regular expression that will match the parts of the filename
        # that come after the sample name
        samp_re = re.compile(r'(_S[0-9]+)?(_L00[1-8])?(_R(1|2))?_001\.((fq(\.gz)?$)|(fastq(\.gz)?$))')
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

    def compile(self, od, wd):
        """Iterate through the dictionary of samples, and fill-in the values
        that were supplied by the user."""
        # For each sample...
        for s in self.samples:
            # Key in the values based on the column names
            self.final_sheet[s] = {
                'FastqR1files': self.samples[s]['R1'],
                'FastqR2file': self.samples[s]['R2'],
                'OutputDir': str(od),
                'WorkingDir': str(wd),
                'TRIM': self.useropts['trim'],
                'trimmomaticOpts': self.finalopts['trimmomatic'],
                'Hisat2index': self.useropts['hisat2_idx'],
                'Hisat2Options': self.useropts['hisat2_threads'] +
                                 self.finalopts['hisat2'],
                'Stranded': self.useropts['stranded'],
                'AnnotationGTF': self.useropts['gtf']
                }
        self.sheet_logger.debug(
            'Samplesheet:\n%s',
            pprint.pformat(self.final_sheet))
        return
