#!/usr/bin/env python

import sys

def sequenceSeg(seq,start,end):
    "Returns a sequence including the lyase, 5000 bp before and 5000 bp after"
    before = start - 5000
    after = end + 5000
    
    shortSeq = seq[before:after]

    return before, after, shortSeq

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
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    writeFile = sys.argv[4]

    seq = readSeq(file)

    before,after,shortSeq = sequenceSeg(seq,start,end)
    targetFile = open(writeFile,'w')
    targetFile.write(file)
    targetFile.write("\n")
    targetFile.write("Start: {0}".format(before))
    targetFile.write("\n")
    targetFile.write("End: {0}".format(after))
    targetFile.write("\n")
    targetFile.write(shortSeq)

    targetFile.close()


if __name__ == "__main__":
    main()

