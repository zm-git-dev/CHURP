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
                if len(tmp) > 2:
                    key = tmp[0] + '\t' + tmp[1]
                    val = tmp[2]
                else:
                    key = tmp[0]
                    val = tmp[1]
                cts_mat[s][key] = val
    # Get a list of the genes
    for sample in sorted(cts_mat):
        genes = [k for k, v in cts_mat[sample].items()]
        genes = sorted(genes)
        break
    # And print them out in order
    print('\t'.join([''] + sorted(cts_mat)))
    for g in genes:
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
