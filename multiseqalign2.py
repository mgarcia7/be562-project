#!/usr/bin/env python

import sys

base_idx = {'A': 0, 'G': 1, 'C': 2, 'T': 3 }
PTR_NONE, PTR_GAP1, PTR_GAP2,PTR_GAP3, PTR_GAP12, PTR_GAP13, PTR_GAP23, PTR_BASE = 0, 1, 2, 3, 4, 5, 6, 7


def multiseqalign2DP(seq1, seq2, seq3, subst_matrix, gap_penalty):
    """
    Return the score of the optimal Needdleman-Wunsch alignment for seq1
    and seq2.
    Note: gap_penalty should be positive (it is subtracted)
    """
    F = [[[0 for k in range(len(seq3)+1)] for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    # Initializes F to matrix of zeros based on length of seq1 and seq2
    TB = [[[PTR_NONE for k in range(len(seq3)+1)] for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    # Initializes TB to matrix of PTR_NONE

    # Initialize dynamic programming table for Needleman-Wunsch alignment
    for i in range(1, len(seq1)+1):
        F[i][0][0] = 0 - i*gap_penalty
        TB[i][0][0] = PTR_GAP3  # indicates a gap in seq3
    for j in range(1, len(seq2)+1):
        F[0][j][0] = 0 - j*gap_penalty
        TB[0][j][0] = PTR_GAP2  # indicates a gap in seq2
    for k in range(1,len(seq3)+1):
        F[0][0][k] = 0 - k*gap_penalty
        TB[0][0][k] = PTR_GAP1  # indicates a gap in seq1
        

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            for k in range(1,len(seq3)+1):
            # Determine index of substitution matrix
            val_i = base_idx[seq1[i-1]]
            val_j = base_idx[seq2[j-1]]
            val_k = base_idx[seq3[k-1]]
            # Determine possible scores and select maximum value
            S1 = F[i-1][j-1][k-1] +  subst_matrix[val_i][val_j][val_k]
            S2 = F[i][j-1][k-1] - gap_penalty
            S3 = F[i-1][j][k-1] - gap_penalty
            S4 = F[i-1][j-1][k] - gap_penalty
            S5 = F[i][j][k-1] - 2*gap_penalty
            S6 = F[i][j-1][k] - 2*gap_penalty
            S7 = F[i-1][j][k] - 2*gap_penalty
            
            F[i][j][k] = max(S1,S2,S3,S4,S5,S6,S7)
            if F[i][j][k] == S1:
                TB[i][j][k] = PTR_BASE
            elif F[i][j][k] == S2:
                TB[i][j][k] = PTR_GAP1
            elif F[i][j][k] == S3:
                TB[i][j][k] = PTR_GAP2
            elif F[i][j][k] == S4:
                TB[i][j][k] = PTR_GAP3
            elif F[i][j][k] == S5:
                TB[i][j][k] = PTR_GAP12
            elif F[i][j][k] == S6:
                TB[i][j][k] = PTR_GAP13               
            else:
                TB[i][j] = PTR_GAP23
    
    return F[len(seq1)][len(seq2)][len(seq3)], F, TB


def traceback(seq1, seq2, seq3, TB):
    s1 = ""
    s2 = ""
    s3 = ""

    i = len(seq1)
    j = len(seq2)
    k = len(seq3)
    

    while TB[i][j][k] != PTR_NONE:
        if TB[i][j][k] == PTR_BASE:
            s1 = seq1[i-1] + s1
            s2 = seq2[j-1] + s2
            s3 = seq3[k-1] + s3
            i = i - 1
            j = j - 1
        elif TB[i][j][k] == PTR_GAP1:
            s1 = '-' + s1
            s2 = seq2[j-1] + s2
            s3 = seq3[k-1] + s3
            j = j - 1
            k = k - 1
        elif TB[i][j][k] == PTR_GAP2:
            s1 = seq1[i-1] + s1
            s2 = '-' + s2
            s3 = seq3[k-1] + s3
            i = i - 1
            k = k - 1
        elif TB[i][j][k] == PTR_GAP3:
            s1 = seq1[i-1] + s1
            s2 = seq2[j-1] + s2
            s3 = '-' + s3
            i = i - 1
            j = j - 1
        elif TB[i][j][k] == PTR_GAP12:
            s1 = '-' + s1
            s2 = '-' + s2
            s3 = seq3[k-1] + s3
            k = k - 1
        elif TB[i][j][k] == PTR_GAP13:
            s1 = '-' + s1
            s2 = seq2[j-1] + s2
            s3 = '-' + s3
            j = j - 1
        elif TB[i][j][k] == PTR_GAP23:
            s1 = seq1[i-1] + s1
            s2 = '-' + s2
            s3 = '-' + s3
            i = i - 1
        else:
            assert False
   
    return s1, s2, s3


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
    [[0, 2, 2, 1], [2, 2, 3, 3], [2, 3, 2, 3], [1, 3, 3, 1]],  
    [2, 2, 3, 3], [2, 0, 1, 2], [2, 1, 1, 3], [3, 2, 3, 2]],    
    [2, 3, 2, 3], [3, 1, 1, 3], [2, 1, 0, 2], [3, 3, 2, 2]],  
    [1, 3, 3, 1], [3, 2, 3, 2], [3, 3, 2, 2], [1, 2, 2, 0]], 
]
gap_penalty = 4


def main():
    # parse command line
    if len(sys.argv) < 3:
        print("Usage: {0} <FASTA 1> <FASTA 2>".format(sys.argv[0]))
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    file3 = sys.argv[3]
    
    #seq1 = sys.argv[1]
    #seq2 = sys.argv[2]
    seq1 = readSeq(file1)
    seq2 = readSeq(file2)
    seq3 = readSeq(file3)

    score, F, TB = multiseqalign2DP(seq1, seq2, seq3, S, gap_penalty)

    print("Score: {0}".format(score))

    s1, s2, s3  = traceback(seq1, seq2, seq3, TB)
    print(s1)
    print(s2)
    print(s3)


if __name__ == "__main__":
    main()
