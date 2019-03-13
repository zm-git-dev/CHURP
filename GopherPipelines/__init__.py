#!/usr/bin/env python
"""GopherPipelines package for building samplesheets and qsub commands to
run analysis on the Minnesota Supercomputing Institute (MSI) cluster."""

# set dunder (double-underscore) variables
__version__ = '0.0.1'

import datetime
import getpass

# Define today's date and the current time
curr_time = datetime.datetime.now()
TODAY = datetime.date.isoformat(datetime.datetime.today())
TIMESTAMP = str(curr_time.hour) + str(curr_time.minute) + str(curr_time.second)
NOW = curr_time.isoformat(sep=' ', timespec='seconds')
UNAME = getpass.getuser()
