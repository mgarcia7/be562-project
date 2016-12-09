import sys
import os
import random
import msa_tools as msa
import genepredictor as gp
import time
import pickle

def main():

	files = ["data/short-sequence/" + n for n in os.listdir("data/short-sequence/")]

	conserved = []

	files = ["data/short-sequence/DshibaeA.txt", "data/short-sequence/JannaschiaCCS1A.txt"]

	for n in range(len(files)-1):
		conserved.append([a for a in msa.Alignment(files[n], files[n+1]).conserved])

	for iterations in range(len(files)):

		a = random.randint(0, len(files)-1)
		b = random.randint(0, len(files)-1)
		while (abs(a - b) < 2):
			b = random.randint(0, len(files)-1)

		conserved.append([a for a in msa.Alignment(files[a], files[b]).conserved])

	#filter results to those of significant length
	conserved = [c for c in conserved if (c.region1.ending - c.region1.beginning > 500)]
	with open('entry.pickle', 'wb') as f:
		for c in conserved:
			pickle.dump(c, f)

	pHMM_dict = gp.initialize_pHMM_models()
	all_score_dicts = gp.score_all_conserved_regions(conserved_regions,pHMM_dicts)
	pred = gp.predictions(conserved_regions,all_score_dicts)


if __name__ == "__main__":
	start_time = time.time()
	main()
	print("--- %s seconds ---" % (time.time() - start_time))