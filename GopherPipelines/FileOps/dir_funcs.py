#!/usr/bin/env python
"""Various simple operations on directories. These will be used in argument
validation functions to make sure we can write into directories required by
the analysis pipeline."""

import pathlib
import tempfile


def dir_exists(d, l):
    """Return true if a directory exists, and false if not."""
    p = pathlib.Path(d)
    l.debug('Checking if %s exists', p)
    if not p.is_dir():
        l.debug('%s does not exist', p)
        return False
    else:
        l.debug('%s exists', p)
        return True


def dir_writeable(d, l):
    """Return true if a directory can be written into, and false if not. We
    use the tempfile library for this: try to actually write a file into the
    directory and then clean it up. We will catch errors and return False if it
    fails."""
    l.debug('Checking if %s is writeable', d)
    try:
        t = tempfile.TemporaryFile(dir=str(d))
        t.close()
        l.debug('%s is writeable', d)
        return True
    except (OSError, PermissionError) as e:
        l.debug('%s is not writeable', d)
        return False


def dir_empty(d, l):
    """Return true if a directory is empty and false if it has files or other
    diretories under it. Note that this will return True if a directory does
    not exist, so it should be verified to exist before being queried here."""
    p = pathlib.Path(d)
    l.debug('Checking if %s is empty', p)
    # get all files and directories under the supplied path. This thankfully is
    # a generator object, so we won't have to worry about the method listing
    # every file on the filesystem if a user tries to supply '/' as the
    # directory
    for f in p.rglob('*'):
        if f:
            l.debug('Found file %s in %s. Not empty', f, p)
            return False
        else:
            continue
    # If we get here, then the generator was empty, so it is a cheap operation
    # anyway
    l.debug('%s appears empty', p)
    return True

def make_dir(d, l):
    """Try to safely make a directory."""
    # Set it to be a pathlib.Path object
    p = pathlib.Path(d)
    l.debug('Making directory %s', p)
    # This should be analogous to mkdir -p
    try:
        p.mkdir(parents=True, exist_ok=True)
        l.debug('Successfully made directory %s', p)
        return True
    except PermissionError as e:
        l.debug('Error making directory %s, permission denied', p)
        return False
    except OSError as e:
        l.debug('Error making directory %s, disk full?', p)
        return False
