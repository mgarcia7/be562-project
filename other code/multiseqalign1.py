#!/usr/bin/env python

import sys

base_idx = {'A': 0, 'G': 1, 'C': 2, 'T': 3 }
PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3


def multiseqalign1DP(seq1, seq2, subst_matrix, gap_penalty):
    """
    Return the score of the optimal Needdleman-Wunsch alignment for seq1
    and seq2.
    Note: gap_penalty should be positive (it is subtracted)
    """
    F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    # Initializes F to matrix of zeros based on length of seq1 and seq2
    TB = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    # Initializes TB to matrix of PTR_NONE

    # initialize dynamic programming table for Needleman-Wunsch alignment
    # (Durbin p.20)
    for i in range(1, len(seq1)+1):
        F[i][0] = 0 - i*gap_penalty
        TB[i][0] = PTR_GAP2  # indicates a gap in seq2
    for j in range(1, len(seq2)+1):
        F[0][j] = 0 - j*gap_penalty
        TB[0][j] = PTR_GAP1  # indicates a gap in seq1

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            # Determine index of substitution matrix
            val_i = base_idx[seq1[i-1]]
            val_j= base_idx[seq2[j-1]]
            # Determine possible scores and select maximum value
            S1 = F[i-1][j-1] +  subst_matrix[val_i][val_j]
            S2 = F[i-1][j] - gap_penalty
            S3 = F[i][j-1] - gap_penalty
            F[i][j] = max(S1,S2,S3)
            if F[i][j] == S1:
                TB[i][j] = PTR_BASE
            elif F[i][j] == S2:
                TB[i][j] = PTR_GAP2
            else:
                TB[i][j] = PTR_GAP1
    
    return F[len(seq1)][len(seq2)], F, TB


def traceback(seq1, seq2, TB):
    s1 = ""
    s2 = ""

    i = len(seq1)
    j = len(seq2)

    while TB[i][j] != PTR_NONE:
        if TB[i][j] == PTR_BASE:
            s1 = seq1[i-1] + s1
            s2 = seq2[j-1] + s2
            i = i - 1
            j = j - 1
        elif TB[i][j] == PTR_GAP1:
            s1 = '-' + s1
            s2 = seq2[j-1] + s2
            j = j - 1
        elif TB[i][j] == PTR_GAP2:
            s1 = seq1[i-1] + s1
            s2 = '-' + s2
            i = i - 1
        else:
            assert False
   
    return s1, s2


def readSeq(filename):
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
    [0, 2, 2, 1],  # A
    [2, 0, 1, 2],  # G
    [2, 1, 0, 2],  # C
    [1, 2, 2, 0]   # T
]
gap_penalty = 4


def main():
    # parse command line
    if len(sys.argv) < 3:
        print("Usage: {0} <FASTA 1> <FASTA 2>".format(sys.argv[0]))
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    
    #seq1 = sys.argv[1]
    #seq2 = sys.argv[2]
    seq1 = readSeq(file1)
    seq2 = readSeq(file2)

    score, F, TB = multiseqalign1DP(seq1, seq2, S, gap_penalty)

    print("Score: {0}".format(score))

    s1, s2  = traceback(seq1, seq2, TB)
    print(s1)
    print(s2)


if __name__ == "__main__":
    main()
