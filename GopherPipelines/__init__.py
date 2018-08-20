#!/usr/bin/env python
"""GopherPipelines package for building samplesheets and qsub commands to
run analysis on the Minnesota Supercomputing Institute (MSI) cluster."""

# set dunder (double-underscore) variables
__version__ = '0.0'

import datetime
import getpass

# Define today's date and the current time
TODAY = datetime.date.isoformat(datetime.datetime.today())
NOW = datetime.datetime.now().isoformat(sep=' ', timespec='seconds')
UNAME = getpass.getuser()

# Define the versions of the software that we are loading on the cluster
# These are the latest available as of 2018-08-20
SW_VERS = {
    'HISAT2': '2.1.0',
    'FASTQC': '0.11.7',
    'TRIMMOMATIC': '0.33',
    'JAVA': 'jdk1.8.0_144',
    'SAMTOOLS': '1.7',
    'R': '3.5.0'
}
