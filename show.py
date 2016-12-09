import sys
import os
import msa_tools as msa

def main():
	files = ["data/short-sequence/DshibaeA.txt", "data/short-sequence/JannaschiaCCS1A.txt"]

	conserved = msa.Alignment(files[0], files[1]).conserved

	print("Alignment  of DShibae and JannaschiaCCS1A:")
	print(conserved[0].region1.organism)
	print(conserved[0].region1.seq)
	print(conserved[0].region2.organism)
	print(conserved[0].region2.seq)

if __name__ == "__main__":
	main()