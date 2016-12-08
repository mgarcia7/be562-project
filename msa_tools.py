#!/usr/bin/env python

import sys
import numpy as np


def compare_conserved(region1, region2):
	base_idx = {'A': 0, 'G': 1, 'C': 2, 'T': 3 }
    PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3
    S = [
        # A   G   C   T
        [ 3, -2, -2, -2],  # A
        [-2,  3, -2, -2],  # G
        [-2, -2,  3, -2],  # C
        [-2, -2, -2,  3] ] # T
    gap_penalty = 4

    F = np.zeros(len(region1.seq) + 1, len(region2.seq), dtype=uint16)
    for i in np.arange(1, F.shape[0]):
    	F[i][0] = F[i-1][0] - gap_penalty
    for i in np.arange(1, F.shape[1]):
    	F[0][i ] F[0][i-1] - gap_penalty

    for i in np.arange(1, F.shape[0]):
    	for j in np.arange(1, F[a].shape[0]):
    		scores = [0, 0, 0]
            scores[0] = F[a-1][b-1] + S[base_idx[seq1[a-1]]][base_idx[seq2[b-1]]]
            scores[1] = F[a][b-1] - gap_penalty
            scores[2] = F[a-1][b] - gap_penalty

            F[a][b] = scores.index(max(scores))

    return F[F.shape[0]-1][F.shape[1]-1]


class Region(object):
    def __init__(self, organism, seq, beginning, ending):
        self.organism = organism
        self.seq = seq
        self.beginning = beginning
        self.ending = ending


class Conserved(object):
    def __init__(self, region1, region2):
        self.region1 = region1
        self.region2 = region2


class Alignment(object):

    def __init__(self, file_1, file_2):
        #assuming 5000bp were allowed on either side of organism's DmdA
        self.organism1, self.seq1 = self.read_fasta(file_1)
        self.organism2, self.seq2 = self.read_fasta(file_2)
        self.seq1_left = self.seq1[:5000]
        self.seq1_right = self.seq1[-5000:]
        self.seq2_left = self.seq2[:5000]
        self.seq2_right = self.seq2[-5000:]
        self.conserved = self.get_conserved(self.organism1, self.organism2,
                                            self.seq1_left, self.seq1_right,
                                            self.seq2_left, self.seq2_right)

    base_idx = {'A': 0, 'G': 1, 'C': 2, 'T': 3 }
    PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3
    S = [
        # A   G   C   T
        [ 3, -2, -2, -2],  # A
        [-2,  3, -2, -2],  # G
        [-2, -2,  3, -2],  # C
        [-2, -2, -2,  3] ] # T
    gap_penalty = 4


    def Smith_Waterman(self, seq1, seq2):
        """
        Return the score of the optimal Smith-Waterman alignment for seq1
        and seq2.
        """
        F = np.zeros([len(seq1) + 1, len(seq2) + 1], dtype=np.uint16)
        TB = np.zeros([len(seq1) + 1, len(seq2) + 1], dtype=np.uint16)

        for a in np.arange(1, F.shape[0]):
            for b in np.arange(1, F[a].shape[0]):
                scores = [0, 0, 0]
                scores[0] = F[a-1][b-1] + self.S[self.base_idx[seq1[a-1]]][self.base_idx[seq2[b-1]]]
                scores[1] = F[a][b-1] - self.gap_penalty
                scores[2] = F[a-1][b] - self.gap_penalty

                maximum = max(scores)
                max_index = scores.index(maximum)
                if max_index == 0:
                    TB[a][b] = self.PTR_BASE
                elif max_index == 1:
                    TB[a][b] = self.PTR_GAP1
                else:
                    TB[a][b] = self.PTR_GAP2

                if maximum < 0:
                    F[a][b] = 0
                else:
                    F[a][b] = maximum

        return F, TB


    def traceback(self, organism1, seq1, d1, organism2, seq2, d2, F, TB):
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
            max_1 = i
            max_2 = j
            while F[i][j] != 0:
                if TB[i][j] == self.PTR_BASE:
                    F[i][j] = 0
                    s1 = seq1[i-1] + s1
                    s2 = seq2[j-1] + s2
                    i = i - 1
                    j = j - 1
                elif TB[i][j] == self.PTR_GAP1:
                    F[i][j] = 0
                    s1 = '-' + s1
                    s2 = seq2[j-1] + s2
                    j = j - 1
                elif TB[i][j] == self.PTR_GAP2:
                    F[i][j] = 0
                    s1 = seq1[i-1] + s1
                    s2 = '-' + s2
                    i = i - 1
                else:
                    assert False

            min_1 = i
            min_2 = j

            #Construct Region objects with name of organism, conserved
            #sequence, and starting and endoing points of that conserved
            #region in that 5000bp flanking region
            #Construct Conserved object with these Regions and append
            region_1 = Region(organism1, s1, min_1 + d1, max_1 + d1)
            region_2 = Region(organism2, s2, min_2 + d2, max_2 + d2)
            conserved.append(Conserved(region_1, region_2))

            #Reset strings used to hold conserved regions and reset
            #new maximum indices to start next traceback
            s1 = ""
            s2 = ""
            m_idx = np.unravel_index(np.argmax(F), F.shape)
            i = m_idx[0]
            j = m_idx[1]

        return conserved


    def get_conserved(self, organism1, organism2, seq1_left, seq1_right, seq2_left, seq2_right):
        #distances passed to traceback are from the left side of the fragment that was passed
        #calculated as such: sequence has length l with organism's copy of lyase
        #with 5000bp on both sides. l - 5000 = length of DmdA. so beginning of conserved
        #regions coming from right flank start at l - 5000 + 5000, or l - 5000 from the
        #left side of fragment.
        conserved = []
        F, TB = self.Smith_Waterman(seq1_left, seq2_left)
        for i in self.traceback(organism1, seq1_left, 0, organism2, seq2_left, 0, F, TB):
            conserved.append(i)
        F, TB = self.Smith_Waterman(seq1_left, seq2_right)
        for i in self.traceback(organism1, seq1_left, 0, organism2, seq2_right, len(self.seq2)-5000, F, TB):
            conserved.append(i)
        F, TB = self.Smith_Waterman(seq1_right, seq2_left)
        for i in self.traceback(organism1, seq1_right, len(self.seq1)-5000, organism2, seq2_left, 0, F, TB):
            conserved.append(i)
        F, TB = self.Smith_Waterman(seq1_right, seq2_right)
        for i in self.traceback(organism1, seq1_right, len(self.seq1)-5000, organism2, seq2_right, len(self.seq2)-5000, F, TB):
            conserved.append(i)
        return conserved


    def read_fasta(self, filename):
        """Reads in a FASTA sequence. Assumes one sequence in the file"""
        seq = []

        with open(filename, "r") as f:
            for line in f:
                if line.startswith(">"):
                    continue
                seq.append(line.rstrip().upper())

        return filename[:(len(filename)-len(".fasta"))], "".join(seq)


"""
def main():
    # parse command line
    if len(sys.argv) < 3:
        print("Usage: {0} <FASTA 1> <FASTA 2>".format(sys.argv[0]))
        sys.exit(1)

    files = [n for n in sys.argv[1:]]

    A = Alignment(sys.argv[1], sys.argv[2])

    print("Alignment:")
    for l in A.conserved:
        print(l.region1.seq)
        print(l.region2.seq)
        print()

    print("Alignment:")
    for l in Alignment(sys.argv[3], sys.argv[4]).conserved:
        print(l.region1.seq)
        print(l.region2.seq)
        print()


if __name__ == "__main__":
    main()
"""
