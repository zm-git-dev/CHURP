#!/usr/bin/env python
"""Various simple operations on directories. These will be used in argument
validation functions to make sure we can write into directories required by
the analysis pipeline."""

import pathlib
import tempfile


def dir_exists(d):
    """Return true if a directory exists, and false if not."""
    p = pathlib.Path(d)
    if not p.is_dir():
        return False
    else:
        return True


def dir_writeable(d):
    """Return true if a directory can be written into, and false if not. We
    use the tempfile library for this: try to actually write a file into the
    directory and then clean it up. We will catch errors and return False if it
    fails."""
    try:
        t = tempfile.TemporaryFile(dir=str(d))
        t.close()
        return True
    except (OSError, PermissionError) as e:
        return False


def dir_empty(d):
    """Return true if a directory is empty and false if it has files or other
    diretories under it."""
    p = pathlib.Path(d)
    # get all files and directories under the supplied path. This thankfully is
    # a generator object, so we won't have to worry about the method listing
    # every file on the filesystem if a user tries to supply '/' as the
    # directory
    for f in p.rglob('*'):
        if f:
            return False
        else:
            continue
    # If we get here, then the generator was empty, so it is a cheap operation
    # anyway
    return True
