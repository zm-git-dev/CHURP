#!/usr/bin/env python
"""Set the default directories based on a base path, the user name, and the
date."""

import pathlib
import getpass
import datetime

# Set the scratch directory as a constant
SCRATCH = pathlib.Path('/panfs/roc/scratch')


def default_output(pipeline):
    """Return a string to a default output scratch directory. This will be
    $SCRATCH/username/YYYY-MM-DD.pipeline_id."""
    # This is ugly, but it is the most flexible way to get the current ISO
    # format date.
    ymd = datetime.date.isoformat(datetime.datetime.today())
    # This is potentially wrong, if the user is running inside a `su' env, or
    # if the user changes some of their environment variables. For the most
    # part, it should be okay.
    uname = getpass.getuser()
    outbase = ymd+'.'+pipeline
    # Path objects are joined with the / operator
    d = SCRATCH / uname / outbase
    return d
