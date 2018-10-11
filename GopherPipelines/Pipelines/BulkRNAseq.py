#!/usr/bin/env python
"""Define a sub-class of the Pipeline class for bulk RNAseq analysis."""

import pprint
import os
import glob
import subprocess
import re
import datetime
import getpass

import GopherPipelines
from GopherPipelines import DieGracefully
from GopherPipelines.Pipelines import Pipeline
from GopherPipelines.SampleSheet import BulkRNASeqSampleSheet
from GopherPipelines.ArgHandling import set_verbosity
from GopherPipelines.FileOps import species_list
from GopherPipelines.FileOps import default_files
from GopherPipelines.FileOps import dir_funcs


class BulkRNAseqPipeline(Pipeline.Pipeline):
    """Sub-class of Pipeline for bulk RNAseq analysis."""

    # Define the pipeline name here
    pipe_name = 'bulk_rnaseq'

    def setup(self, args):
        """Initialize the pipeline object. We will call the general
        Pipeline.__init__() here, as well as set some specific pipeline
        attributes."""
        # Set up the verbosity and logging
        self.pipe_logger = set_verbosity.verb(args['verbosity'], __name__)
        # First validate the arguments. We do this in the subclass because the
        # dependencies of the arguments are pipeline-specific.
        valid_args = self._validate_args(args)
        self.pipe_logger.debug('New BulkRNAseqPipeline instance.')
        self.pipe_logger.debug(
            'Validated args:\n%s', pprint.pformat(valid_args))

        # Set the minimum gene length
        self.min_gene_len = str(valid_args['mingene'])

        # And make a sample sheet from the args
        self.sheet = BulkRNASeqSampleSheet.BulkRNASeqSampleSheet(valid_args)

        # Set the paths to the single sample analysis script. This will be
        # submitted to the scheduler. This is a little ugly, but because the
        # package contains the shell script data, it should be OK to define a
        # relative path. __file__ is the currently running script,
        # os.path.realpath() gives the full path, with symlinks resolved. We
        # then rsplit() the string returned by realpath() to get the base dir
        # of the gopher-pipelines scipt. Woof.
        self.single_sample_script = os.path.join(
            os.path.realpath(__file__).rsplit(os.path.sep, 3)[0],
            'PBS',
            'bulk_rnaseq_single_sample.pbs')
        self.summary_script = os.path.join(
            os.path.realpath(__file__).rsplit(os.path.sep, 3)[0],
            'PBS',
            'run_summary_stats.pbs')
        self.de_script = os.path.join(
            os.path.realpath(__file__).rsplit(os.path.sep, 3)[0],
            'R_Scripts',
            'summarize_bulk_rnaseq.R')
        return

    def _validate_args(self, a):
        """Validate arguments for the BulkRNAseqPipeline object. We define it
        in this file because it only really needs to be accessible to this
        subclass. Argument dependencies are pipeline-specific. For the bulk
        RNAseq analysis pipeline:
            - FASTQ dir and UMGC sheet are mutually exclusive
            - Organism and HISAT2 index + GTF are mutually exclusive
            - Check that the HISAT2 index is complete
        Further, sanitize the paths of the output dir, working dir, and hisat2
        index.
        """
        # Check the completeness of the argument dictionary. -f must be
        # specified and (-x and -g) or -r must be specified. After checking the
        # FASTQ folder, check the helper commands
        if not a['fq_folder']:
            DieGracefully.die_gracefully(DieGracefully.BRNASEQ_INC_ARGS)
        elif not ((a['hisat2_idx'] and a['gtf']) or a['organism']):
            DieGracefully.die_gracefully(DieGracefully.BRNASEQ_INC_ARGS)
        elif (a['hisat2_idx'] and a['gtf']) and a['organism']:
            DieGracefully.die_gracefully(DieGracefully.BRNASEQ_CONFLICT)
        # Convert all of the paths into absolute paths
        a['fq_folder'] = os.path.realpath(
            os.path.expanduser(str(a['fq_folder'])))
        a['hisat2_idx'] = os.path.realpath(
            os.path.expanduser(str(a['hisat2_idx'])))
        a['gtf'] = os.path.realpath(os.path.expanduser(str(a['gtf'])))
        a['outdir'] = os.path.realpath(os.path.expanduser(str(a['outdir'])))
        a['workdir'] = os.path.realpath(os.path.expanduser(str(a['workdir'])))
        if a['expr_groups']:
            a['expr_groups'] = os.path.realpath(os.path.expanduser(str(
                a['expr_groups'])))
        self.pipe_logger.debug('GTF: %s', a['gtf'])
        self.pipe_logger.debug('Adapters: %s', a['adapters'])
        self.pipe_logger.debug('FASTQ Folder: %s', a['fq_folder'])
        self.pipe_logger.debug('Output Dir: %s', a['outdir'])
        self.pipe_logger.debug('Working Dir: %s', a['workdir'])
        self.pipe_logger.debug('HISAT2 Idx: %s', a['hisat2_idx'])
        self.pipe_logger.debug('Expr Groups: %s', a['expr_groups'])
        # Check that the adapters and GTF file exist
        try:
            handle = open(a['gtf'], 'r')
            handle.close()
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.BAD_GTF)
        if not a['adapters']:
            a['adapters'] = '$TRIMMOMATIC/adapters/all_illumina_adapters.fa'
        else:
            try:
                a['adapters'] = os.path.realpath(
                    os.path.expanduser(str(a['adapters'])))
                handle = open(a['adapters'], 'r')
                handle.close()
            except OSError:
                DieGracefully.die_gracefully(DieGracefully.BAD_ADAPT)
        if a['expr_groups']:
            try:
                handle = open(a['expr_groups'], 'r')
                handle.close()
            except OSError:
                DieGracefully.die_gracefully(DieGracefully.BRNASEQ_BAD_GPS)
        # Validate the FASTQ folder
        self._validate_fastq_folder(a['fq_folder'])
        # Validate the hisat2 index
        self._validate_hisat_idx(a['hisat2_idx'])
        return a

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

    def _validate_hisat_idx(self, i):
        """Raise an error if the provided HISAT2 index is not complete -
        all of the [1-8].ht2l? files should be present."""
        # Build glob patterns for the normal and long indices
        norm_idx = i + '.[1-8].ht2'
        long_idx = i + '.[1-8].ht2l'
        # Do the search
        self.pipe_logger.debug('Searching for %s', norm_idx)
        norm_idx_files = glob.glob(norm_idx)
        self.pipe_logger.debug('Found %i idx files', len(norm_idx_files))
        # There should be 8 total
        if len(norm_idx_files) == 8:
            return
        else:
            self.pipe_logger.debug(
                'Normal idx not found. Searching for long idx.')
            long_idx_files = glob.glob(long_idx)
            self.pipe_logger.debug(
                'Found %i long idx files', len(long_idx_files))
            if len(long_idx_files) == 8:
                return
            else:
                self.pipe_logger.error('Cound not find HISAT2 idx files!')
                DieGracefully.die_gracefully(DieGracefully.BAD_HISAT)
        return

    def _prepare_samplesheet(self):
        """Call the samplesheet build method here. This will build the
        dictionary that will hold all samplesheet data, and then write it into
        the output directory."""
        self._run_checks()
        self.sheet.compile(self.outdir, self.workdir)
        ss_path = self.sheet.write_sheet(self.outdir, self.pipe_name, '|')
        return ss_path

    def qsub(self):
        """Write the qsub command. We will need the path to the samplesheet,
        the number of samples in the samplesheet, and the scheduler options
        that were passed as arguments. This is defined in the subclass rather
        than in the main Pipeline class because the exact form of the qsub
        command depends on which pipeline we are running."""
        ss = self._prepare_samplesheet()
        # Make the script filename
        pname = default_files.default_pipeline(self.pipe_name)
        pname = os.path.join(self.outdir, pname)
        if os.path.isfile(pname):
            self.pipe_logger.warning(
                'Submission script %s already exists. Overwriting!', pname)
        try:
            handle = open(pname, 'w')
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.BAD_OUTDIR)
        # Write the header of the script
        handle.write('#!/bin/bash\n')
        # Write some meta-data lines
        handle.write('# Generated by gopher-pipelines version ' +
                     GopherPipelines.__version__ + '\n')
        handle.write('# Generated at ' + GopherPipelines.NOW + '\n')
        self.pipe_logger.debug(
            'Number of samples: %i', len(self.sheet.final_sheet))
        self.pipe_logger.debug('Samplesheet: %s', ss)
        # Write command to figure out email address of submitting user
        handle.write('user_name="$(id -u -n)"\n')
        handle.write('user_email="${user_name}@umn.edu"\n')
        # This is the command for aligning and cleaning
        qsub_resources = '"mem=' + str(self.mem) + 'mb'
        qsub_resources += ',nodes=1:ppn=' + str(self.ppn)
        qsub_resources += ',walltime=' + str(self.walltime * 3600) + '"'
        # Set the group string here
        if self.group:
            qsub_group = '-A' + ' ' + self.group
        else:
            qsub_group = ''
        qsub_array = '1'
        if len(self.sheet.final_sheet) > 1:
            qsub_array += '-' + str(len(self.sheet.final_sheet))
        # Write a few variables into the header of the script so they are
        # easy to find
        handle.write('OUTDIR=' + '"' + self.outdir + '"\n')
        handle.write('WORKDIR=' + '"' + self.workdir + '"\n')
        handle.write('DE_SCRIPT=' + '"' + self.de_script + '"\n')
        handle.write('SAMPLESHEET=' + '"' + ss + '"\n')
        handle.write('PURGE=' + '"' + self.purge + '"\n')
        aln_cmd = [
            'qsub',
            '-q', 'mesabi',
            '-m', 'abe',
            '-M', '"${user_email}"',
            qsub_group,
            '-o', '"${OUTDIR}"',
            '-e', '"${OUTDIR}"',
            '-l', qsub_resources,
            '-t', qsub_array,
            '-v', '"SampleSheet=${SAMPLESHEET},PURGE=${PURGE}"',
            self.single_sample_script]
        # Write the first qsub command
        handle.write('single_id=$(' + ' '.join(aln_cmd) + ')\n')
        # This is the command for counting and normalizing reads
        summary_vars = ''.join([
            'SampleSheet=${SAMPLESHEET}',
            ',MINLEN=',
            self.min_gene_len,
            ',RSUMMARY=${DE_SCRIPT}'])
        summary_cmd = [
            'qsub',
            '-q', 'mesabi',
            '-m', 'abe',
            '-M', '"${user_email}"',
            qsub_group,
            '-o', '"${OUTDIR}"',
            '-e', '"${OUTDIR}"',
            '-W', '"depend=afterokarray:${single_id}"',
            '-v', '"' + summary_vars + '"',
            self.summary_script]
        # Write the second command
        handle.write('summary_id=$(' + ' '.join(summary_cmd) + ')\n')
        # Write some echo statements for users' information
        handle.write('echo "Output and logs will be written to ${OUTDIR}"\n')
        handle.write('echo "Emails will be sent to ${user_email}"\n')
        handle.write('echo "Single samples job array ID: ${single_id}"\n')
        handle.write('echo "Summary job ID: ${summary_id}"\n')
        self.pipe_logger.debug('qsub:\n%s', ' '.join(aln_cmd))
        self.pipe_logger.debug('qsub:\n%s', ' '.join(summary_cmd))
        handle.flush()
        handle.close()
        return (pname, ss)
