#!/usr/bin/env python
"""Script to merge a bunch of individual counts files from HTSeq. The final
matrix will have samples as columns and genes as rows. Takes one argument:
    1) Directory with HTSeq counts files
"""

import sys
import os
import glob


def find_counts_files(c):
    """Return a list of the counts files found in the supplied directory."""
    cfiles = {}
    apath = os.path.abspath(os.path.expanduser(c))
    for cf in glob.iglob(apath + '/**/*.counts'):
        sn = os.path.basename(cf).rstrip('.counts')
        cfiles[sn] = cf
    return cfiles


def main(counts_dir):
    """Main funcction."""
    cts_files = find_counts_files(counts_dir)
    # Then, store the counts in a dictionary, keyed on sample name
    cts_mat = {}
    for s in cts_files:
        cts_mat[s] = {}
        with open(cts_files[s], 'r') as f:
            for line in f:
                tmp = line.strip().split()
                cts_mat[s][(tmp[0], tmp[1])] = tmp[2]
    # And print them out in order
    print('\t'.join([''] + sorted(cts_mat)))
    for g in sorted(cts_mat[sorted(cts_mat)[0]]):
        toprint = [g]
        for s in sorted(cts_mat):
            toprint.append(cts_mat[s][g])
        print('\t'.join(toprint))
    return


if len(sys.argv) != 2:
    print("""Script to merge a bunch of individual counts files from HTSeq. The final
matrix will have samples as columns and genes as rows. Takes one argument:
    1) Directory with HTSeq counts files""")
    exit(1)
else:
    main(sys.argv[1])
