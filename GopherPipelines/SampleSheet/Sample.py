#!/usr/bin/env python
"""Define a class for Sample handling/creation."""

class Samplesheet(object):
	"""We set the following attributes:
		self.name
		self.r1
		self.r2
		self.tech
		self.trim_opts
		self.hisat2_opts
	
	And define the following methods:
		resolve_options(): Handles user-provided options
		sanitize_paths(): Set full paths for I/O.
	"""
	
	def __init__(self, args):
		"""Initialize the object with relevant info that is parsed elsewhere like:
			"""
		self.name = []
		self.r1 = []
		self.r2 = []
		self.tech = []
		self.trim_opts = []
		self.hisat2_opts = []
	
	def resolve_options(self):
		""" Given some set of options, do some things that need doing
	"""
		pass
		
		
	def santize_paths(self):
		"""Make paths sane, create absolute paths"
	"""
		pass