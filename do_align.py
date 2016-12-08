import sys
import random
import msa_tools

def main():

	# parse command line
    if len(sys.argv) < 3:
        print("Usage: {0} <FASTA 1> <FASTA 2>".format(sys.argv[0]))
        sys.exit(1)

    files = [n for n in sys.argv[1:]]

    conserved = []

    for n in range(len(files)-1):
    	conserved.append(Alignment(files[n], files[n+1]).conserved)

    for iterations in range(len(files)):

	    a = random.randint(0, len(files)-1)
	    b = random.randint(0, len(files)-1)
	    while (abs(a - b) < 2):
	    	b = random.randint(0, len(files)-1)

	    conserved.append(Alignment(files[a], files[b]).conserved)



if __name__ == "__main__":
	main()