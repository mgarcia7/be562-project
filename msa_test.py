#!/usr/bin/env python

import sys
import numpy as np

base_idx = {'A': 0, 'G': 1, 'C': 2, 'T': 3 }
PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3

def seqalignDP(seq1, seq2, subst_matrix, gap_penalty):
    """
    Return the score of the optimal Smith-Waterman alignment for seq1
    and seq2.
    """

    #F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    #TB = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    F = np.zeros([len(seq1) + 1, len(seq2) + 1], dtype=np.uint32)
    TB = np.zeros([len(seq1) + 1, len(seq2) + 1], dtype=np.uint32)

    # YOUR CODE HERE
    # Fill in the dynamic programming tables F and TB, starting at [1][1]
    # Hints: The first row and first column of the table F[i][0] and F[0][j]
    # contain dummy values.
    #  (see for illustration Durbin p.21, Figure 2.5, but be careful what you
    #   think of as rows and what you think of as columns)
    #  Hence, the bases corresponding to F[i][j] are actually seq1[i-1] and
    #  seq2[j-1].
    #  Use the dictionary base_idx to convert from the character to an index to
    #   look up entries of the substitution matrix.
    #  To get started, you can complete and run the algorithm filling in only
    #    F, and then figure out how to do TB.

    for a in np.arange(1, F.shape[0]):
        for b in np.arange(1, F[a].shape[0]):
            scores = [0, 0, 0]
            scores[0] = F[a-1][b-1] + S[base_idx[seq1[a-1]]][base_idx[seq2[b-1]]]
            scores[1] = F[a][b-1] - gap_penalty
            scores[2] = F[a-1][b] - gap_penalty

            maximum = max(scores)
            max_index = scores.index(maximum)
            if max_index == 0:
                TB[a][b] = PTR_BASE
            elif max_index == 1:
                TB[a][b] = PTR_GAP1
            else:
                TB[a][b] = PTR_GAP2

            if maximum < 0:
                F[a][b] = 0
            else:
                F[a][b] = maximum

    return F, TB


def traceback(seq1, seq2, F, TB):
    max_iter = 5
    iterations = 0

    s1 = ""
    s2 = ""

    conserved = []

    m_idx = np.unravel_index(np.argmax(F), F.shape)
    i = m_idx[0]
    j = m_idx[1]

    while iterations < max_iter:
        iterations += 1
        while F[i][j] != 0:
            if TB[i][j] == PTR_BASE:
                F[i][j] = 0
                s1 = seq1[i-1] + s1
                s2 = seq2[j-1] + s2
                i = i - 1
                j = j - 1
            elif TB[i][j] == PTR_GAP1:
                F[i][j] = 0
                s1 = '-' + s1
                s2 = seq2[j-1] + s2
                j = j - 1
            elif TB[i][j] == PTR_GAP2:
                F[i][j] = 0
                s1 = seq1[i-1] + s1
                s2 = '-' + s2
                i = i - 1
            else:
                assert False

        conserved.append(s1)
        conserved.append(s2)
        s1 = ""
        s2 = ""

        m_idx = np.unravel_index(np.argmax(F), F.shape)
        i = m_idx[0]
        j = m_idx[1]

    return conserved

def print_matrix(mat):
    for a in mat:
        print(a)

def read_fasta(filename):
    """Reads in a FASTA sequence. Assumes one sequence in the file"""
    seq = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.rstrip().upper())
    return "".join(seq)


# Substituation matrix and gap_penalty
S = [
    # A   G   C   T
    [ 3, -1, -2, -2],  # A
    [-1,  3, -2, -2],  # G
    [-2, -2,  3, -1],  # C
    [-2, -2, -1,  3]   # T
]
gap_penalty = 4


def main():
    # parse command line
    if len(sys.argv) < 3:
        print("Usage: {0} <FASTA 1> <FASTA 2>".format(sys.argv[0]))
        sys.exit(1)

    files = [n for n in sys.argv[1:]]

    genomes = [read_fasta(n) for n in files]

    F, TB = seqalignDP(genomes[0], genomes[1], S, gap_penalty)

    conserved = []

    conserved.append(traceback(genomes[0], genomes[1], F, TB))

    print("Alignment:")
    for l in conserved:
        for s in l:
            print(s)
            print()
    #print("F:")
    #print_matrix(F)
    #print("TB:")
    #print_matrix(TB)

if __name__ == "__main__":
    main()
