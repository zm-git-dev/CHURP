#!/usr/bin/env python
"""Return a logger object that will be used throughout the pipeline to print
messages to the terminal. The level of logging is determined by user input."""

import logging

def verb(level, lname):
    """Return a logging.getLogger object that will print messages to the
    console for debug purposes."""
    user_levels = {
        'debug': logging.DEBUG,
        'info': logging.INFO,
        'warn': logging.WARNING
        }
    # Start a new logger, and give it the name of the current module
    l = logging.getLogger(lname)
    l.setLevel(user_levels[level])
    # Add a stream handler, so messages will be printed to the console
    s = logging.StreamHandler()
    # Set the log level
    s.setLevel(user_levels[level])
    # Set the log format:
    #   time - module_name - level: message
    fmt = '%(asctime)s - %(name)s - %(levelname)s: %(message)s'
    formatter = logging.Formatter(fmt)
    s.setFormatter(formatter)
    # Add the handler to the logger
    l.addHandler(s)
    # Return the logging object for use
    return l
