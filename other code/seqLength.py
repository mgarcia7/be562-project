#!/usr/bin/env python

import sys

def seqLength(sequence):
    # Returns the length of the input sequence
    seqL = len(sequence)+1

    return seqL

def readSeq(filename):
    """Reads in a FASTA sequence. Assumes one sequence in the file"""
    seq = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.rstrip().upper())

    return "".join(seq)

def main():

    file = sys.argv[1]
    sequence = readSeq(file)
    #seqL = seqLength(sequence)
    seqL = len(sequence)+1
    print(seqL)
    print(sequence)
    #print('Sequence length :{0}'.format(seqL))

    if __name__ == "__main__":
        main()
