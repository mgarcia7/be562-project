import sys
import os
import random
import msa_tools

def main():

	# parse command line
    #if len(sys.argv) < 3:
    #    print("Usage: {0} <FASTA 1> <FASTA 2>".format(sys.argv[0]))
    #    sys.exit(1)

    files = ["/data/genomes/" + n for n in os.listdir("/data/genomes")]

    conserved = []

    for n in range(len(files)-1):
    	conserved.append(Alignment(files[n], files[n+1]).conserved)

    for iterations in range(len(files)):

	    a = random.randint(0, len(files)-1)
	    b = random.randint(0, len(files)-1)
	    while (abs(a - b) < 2):
	    	b = random.randint(0, len(files)-1)

	    conserved.append(Alignment(files[a], files[b]).conserved)

	#filter results to those of significant length
	conserved = [c for c in conserved if (c.seq1.ending - c.seq1.beginning > 500)]

if __name__ == "__main__":
	main()