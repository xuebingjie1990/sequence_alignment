#!/usr/bin/env python

"""
    usage:
        align_sequences [options] seq1.fa seq2.fa
    where the options are:
        -h,--help : print usage and quit
        -m,--match: score of a match in the alignment [2]
        -x,--mismatch: penalty for a mismatch in the alignment [1]
        -g,--gapopen: penalty for opening a new gap [4]
        -e,--gapextend: penalty for extending a gap [1]
"""

from sys import argv, stderr
from getopt import getopt, GetoptError
import numpy as np
from itertools import groupby

# a simple function to read the name and sequence from a file
# The file is expected to have just one contig/sequence. This function
# checks the assumption and complains if it is not the case.
def read_single_contig_fasta(filename):
    names = []
    sequences = []
    with open(filename, "r") as f:
        line = f.readline()
        assert line.startswith(">")
        names.append(line.strip().split("\t"))
        sequence = ""
        for line in f:
            if line.startswith(">"):
                sequences.append(sequence)
                names.append(line.strip().split("\t"))
                sequence = ""
            else:
                for x in line.strip():
                    if x not in ["A", "C", "G", "T"]:
                        print("Unknown nucleotide {}".format(x), file=stderr)
                        exit(3)
                sequence += line.strip()

    sequences.append(sequence)
    assert len(names) == 1
    assert len(sequences) == 1
    return names[0], sequences[0]


def smith_waterman(seq1, seq2, match, mismatch, gapopen, gapextend):
    # Initialising the variables
    max_score = 0
    max_index = (-1, -1)
    alnseq1 = ""
    alnseq2 = ""

    # Generating the empty scoring and tracing matrices
    row = len(seq1) + 1
    col = len(seq2) + 1
    scoring_matrix = np.zeros(shape=(row, col), dtype=np.float64)
    tracing_matrix = np.zeros(shape=(row, col), dtype=np.float64)

    # where the cell's value is coming from
    stop = 0
    left = 1
    up = 2
    diagonal = 3

    # Calculating the scores for all cells in the matrix
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score (match score)
            match_value = match if seq1[i - 1] == seq2[j - 1] else -mismatch
            diagonal_score = scoring_matrix[i - 1, j - 1] + match_value

            # Calculating the vertical gap score
            if tracing_matrix[i - 1, j] == up:
                count_dups = [
                    sum(1 for _ in group)
                    for _, group in groupby(scoring_matrix[0 : i - 1, j])
                ]
                numgap = count_dups[-1]
                vertical_score = (
                    scoring_matrix[i - numgap, j] - (numgap + 1) * gapextend
                )
            else:
                vertical_score = scoring_matrix[i - 1, j] - gapopen

            # Calculating the horizontal gap score
            if tracing_matrix[i, j - 1] == left:
                count_dups = [
                    sum(1 for _ in group)
                    for _, group in groupby(scoring_matrix[i, 0 : j - 1])
                ]
                numgap = count_dups[-1]
                horizontal_score = (
                    scoring_matrix[i, j - numgap] - (numgap + 1) * gapextend
                )
            else:
                horizontal_score = scoring_matrix[i, j - 1] - gapopen

            # Taking the highest score
            scoring_matrix[i, j] = max(
                0, diagonal_score, vertical_score, horizontal_score
            )

            # Tracking where the cell's value is coming from
            if scoring_matrix[i, j] == diagonal_score:
                tracing_matrix[i, j] = diagonal

            elif scoring_matrix[i, j] == horizontal_score:
                tracing_matrix[i, j] = left

            elif scoring_matrix[i, j] == vertical_score:
                tracing_matrix[i, j] = up

            elif scoring_matrix[i, j] == 0:
                tracing_matrix[i, j] = stop

            # Tracking the cell with the maximum score
            if scoring_matrix[i, j] >= max_score:
                max_index = (i, j)
                max_score = scoring_matrix[i, j]

    # Initialising the variables for tracing
    current_alnseq1 = ""
    current_alnseq2 = ""
    (max_i, max_j) = max_index

    # Tracing and computing the pathway with the local alignment
    while tracing_matrix[max_i, max_j] != stop:
        if tracing_matrix[max_i, max_j] == diagonal:
            current_alnseq1 = seq1[max_i - 1]
            current_alnseq2 = seq2[max_j - 1]
            max_i = max_i - 1
            max_j = max_j - 1

        elif tracing_matrix[max_i, max_j] == up:
            current_alnseq1 = seq1[max_i - 1]
            current_alnseq2 = "-"
            max_i = max_i - 1

        elif tracing_matrix[max_i, max_j] == left:
            current_alnseq1 = "-"
            current_alnseq2 = seq2[max_j - 1]
            max_j = max_j - 1

        alnseq1 = alnseq1 + current_alnseq1
        alnseq2 = alnseq2 + current_alnseq2

    # Reversing the order of the sequences
    alnseq1 = alnseq1[::-1]
    alnseq2 = alnseq2[::-1]

    return max_score, alnseq1, alnseq2


def main(filename1, filename2, match, mismatch, gapopen, gapextend):
    # read the name and sequence from the file
    name1, seq1 = read_single_contig_fasta(filename1)
    name2, seq2 = read_single_contig_fasta(filename2)

    # this function takes as input two nucleotide sequences along with
    # scores for an alignment match, mismatch, opening a new gap, and
    # extending an existing gap. This should return the maximum alignment
    # score as well as the alignment. For examples see the testdriver script
    max_score, alnseq1, alnseq2 = smith_waterman(
        seq1, seq2, match, mismatch, gapopen, gapextend
    )

    print("Maximum alignment score: {}".format(max_score))
    print("Sequence1 : {}".format(alnseq1))
    print("Sequence2 : {}".format(alnseq2))


if __name__ == "__main__":
    try:
        opts, args = getopt(
            argv[1:],
            "hm:x:g:e:",
            ["help", "match=", "mismatch=", "gapopen=", "gapextend="],
        )
    except GetoptError as err:
        print(err)
        print(__doc__, file=stderr)
        exit(1)

    match = 2
    mismatch = 1
    gapopen = 4
    gapextend = 1

    for o, a in opts:
        if o in ("-h", "--help"):
            print(__doc__, file=stderr)
            exit()
        elif o in ("-m", "--match"):
            match = float(a)
        elif o in ("-x", "--mismatch"):
            mismatch = float(a)
        elif o in ("-g", "--gapopen"):
            gapopen = float(a)
        elif o in ("-e", "--gapextend"):
            gapextend = float(a)
        else:
            assert False, "unhandled option"

    if len(args) != 2:
        print(__doc__, file=stderr)
        exit(2)

    main(args[0], args[1], match, mismatch, gapopen, gapextend)
