#!/usr/bin/env python
"""Define a class for Sample sheet handling/creation."""

import glob, sys, os
from GopherPipelines.FileOps import default_dirs
from GopherPipelines.ArgHandling import bulk_rnaseq_args
from GopherPipelines.ArgHandling import sc_rnaseq_args

class Samplesheet(object):
	"""We set the following attributes:
		self.blah
		self.blah2
		self.blah3
		sheet_dict = {}
	
	And define the following methods:
		parse_sheet(): Handles a user-provided sample sheet
		create_sheet(): Creates a sample sheet from UMGC data dump.
	"""
	
	def __init__(self, args):
		"""Initialize the object with relevant info that is parsed elsewhere like:
			- path to fastq"""
		self.fq_dir = args['fq_folder']
		self.sheet_dict = {}
	def parse_sheet(ss):
		"""Method for parsing a provided sample sheet, i.e., if ss in Pipelines.prepare_samplesheet()"""
		with open(ss) as sheet:
			for line in sheet:
				add_relevant_stuff_to_dict
	
	def create_sheet(ss):
		"""Method for creating a sample sheet from scratch given the UMGC-provided data"""
		with open(UMGC_spreadsheetloc) as sheet:
			add_stuff_to_dict
