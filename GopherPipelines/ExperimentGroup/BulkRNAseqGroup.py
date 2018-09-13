#!/usr/bin/env python
"""Experiment group for the bulk RNAseq pipeline."""

import pprint
import os

import GopherPipelines
from GopherPipelines import DieGracefully
from GopherPipelines.ExperimentGroup import ExpGroup
from GopherPipelines.ArgHandling import set_verbosity


class BulkRNAseqGroup(ExpGroup.ExpGroup):
    """Inherits from the ExpGroup object. Some of these functions will be very
    similar to those that are written for the bulk RNAseq samplesheet object.
    This is unavoidable, unfortunately."""
    pass

    def setup(self, args):
        """Do a detailed check of the arguments passed to this function. We
        want to validate the output directory, working directory, FASTQ
        directory, and the columns that were passed."""
        valid_args = self._validate(args)
        # Append the extra columns to the default ones
        self.columns.extend(valid_args['extra_column'])
        self.samples = self._get_sample_names(valid_args['fq_folder'])
        self._build_groups()
        return

    def _validate(self, a):
        """Validate the arguments. We want to make sure that the FASTQ
        directory is not empty, the columns do not collide with each other, and
        that the names do not have any commas in them."""
        self._validate_fastq_folder(a['fq_folder'])
        # Drop a warning that specifying extra columns means that there will be
        # some more specialized statistical analysis required
        if a['extra_column']:
            self.group_logger.warning(
                'Specifying additional columns for experimental conditions '
                'is an advanced feature, and will require you to write custom '
                'scripts for statistical analysis. gopher-pipelines will do '
                'tests on the "Group" column (present by default), but will '
                'not account for additional experimental details in your '
                'design. This is not an error message.')
        # Check the experimental columns - first make sure that the names are
        # not duplicated
        tot_col = self.columns + a['extra_column']
        if len(tot_col) != len(set(tot_col)):
            self.group_logger.warning(
                'Duplicate columns specified. This will not cause an error ' +
                'in the Python script, but it may cause an error in any ' +
                'downstream statistical analysis.')
        # Check the supplied columns for bad values
        for e in a['extra_column']:
            if ',' in e:
                self.group_logger.error('Column names cannot contain commas.')
                DieGracefully.die_gracefully(DieGracefully.GROUP_BAD_COL)
        # Turn relative paths into absolute paths
        a['fq_folder'] = os.path.realpath(
            os.path.expanduser(a['fq_folder']))
        return a
