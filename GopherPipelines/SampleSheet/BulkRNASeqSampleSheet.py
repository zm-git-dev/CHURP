#!/usr/bin/env python
"""Define a sub-class of samplesheet that holds the data for bulk RNAseq
samples."""

import os
import re
import pprint

from GopherPipelines import DieGracefully
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
        if args['rmdup']:
            self.useropts['rmdup'] = 'no'
        else:
            self.useropts['rmdup'] = 'yes'
        if args['unstranded']:
            self.useropts['unstranded'] = 'yes'
        else:
            self.useropts['unstranded'] = 'no'
        self.useropts['gtf'] = args['gtf']
        self.useropts['hisat2_idx'] = args['hisat2_idx']
        self.useropts['hisat2_threads'] = '-p ' + str(args['ppn']) + ' '
        # Set the column order to be the columns of the sample sheet. This will
        # eventually become the header of the sheet.
        self.column_order.extend([
            'Group',
            'FastqR1files',
            'FastqR2file',
            'OutputDir',
            'WorkingDir',
            'TRIM',
            'RMDUP',
            'trimmomaticOpts',
            'Hisat2index',
            'Hisat2Options',
            'Unstranded',
            'AnnotationGTF'])
        self._get_fq_paths(args['fq_folder'])
        self._resolve_options()
        # Set the sample group memberships based on the expr_groups argument
        self._set_groups(args['expr_groups'])
        return

    def _set_groups(self, groups):
        """If the groups argument is NoneType, then there were no experimental
        groups passed to the pipeline script, and we fill in 'NULL' for each
        group. If it was passed, then we parse it for group memberships. We set
        any overlapping sample groups to be the same value. We report any
        non-overlapping samples as warnings, and complete non-overlap as an
        error."""
        if not groups:
            self.sheet_logger.debug(
                'No groups file passed. All samples are NULL group.')
            for s in self.samples:
                self.samples[s]['Group'] = 'NULL'
            return
        else:
            self.sheet_logger.debug('Parsing %s for groups.', groups)
            csv_gps = {}
            with open(groups, 'r') as f:
                for index, line in enumerate(f):
                    if index == 0:
                        continue
                    else:
                        tmp = line.strip().split(',')
                        csv_gps[tmp[0]] = tmp[1]
            self.sheet_logger.debug(
                'CSV experimental groups:\n%s',
                pprint.pformat(csv_gps))
            # Calculate the overlaps
            fq_samples = set(self.samples)
            csv_samples = set(csv_gps)
            fq_only = fq_samples - csv_samples
            csv_only = csv_samples - fq_samples
            overlap = fq_samples & csv_samples
            # Thow an error if there are no overlapping samples
            if len(overlap) == 0:
                self.sheet_logger.error(
                    'The FASTQ directory and CSV file appear to mismatch.')
                DieGracefully.die_gracefully(DieGracefully.BRNASEQ_NO_SAMP_GPS)
            # Drop warnings if there are exclusive samples
            if len(fq_only) > 0:
                self.sheet_logger.warning(
                    'Samples found that do not have CSV entries, setting to '
                    'NULL group:\n%s',
                    '\n'.join(list(fq_only)))
            if len(csv_only) > 0:
                self.sheet_logger.warning(
                    'Ignoring samples with groups but no reads in the FASTQ '
                    'directory:\n%s',
                    '\n'.join(list(csv_only)))
            # Iterate through the sample dictionary and set groups
            for s in self.samples:
                gp = csv_gps.get(s, 'NULL')
                self.sheet_logger.debug('Sample %s gets group %s.', s, gp)
                self.samples[s]['Group'] = gp
            return

    def compile(self, od, wd):
        """Iterate through the dictionary of samples, and fill-in the values
        that were supplied by the user."""
        # For each sample...
        for s in self.samples:
            # Key in the values based on the column names
            self.final_sheet[s] = {
                'Group': self.samples[s]['Group'],
                'FastqR1files': self.samples[s]['R1'],
                'FastqR2file': self.samples[s]['R2'],
                'OutputDir': str(od),
                'WorkingDir': str(wd),
                'TRIM': self.useropts['trim'],
                'RMDUP': self.useropts['rmdup'],
                'trimmomaticOpts': self.finalopts['trimmomatic'],
                'Hisat2index': self.useropts['hisat2_idx'],
                'Hisat2Options': self.useropts['hisat2_threads'] +
                                 self.finalopts['hisat2'],
                'Unstranded': self.useropts['unstranded'],
                'AnnotationGTF': self.useropts['gtf']
                }
        self.sheet_logger.debug(
            'Samplesheet:\n%s',
            pprint.pformat(self.final_sheet))
        return
