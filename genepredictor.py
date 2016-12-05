################################################################
# TASK: Predict what gene it is based on annotated genes + BLAST
# INPUT: Conserved regions
# OUTPUT: Name of gene of all the conserved regions/likelihood that it is that gene
################################################################
import numpy as np
import re
from collections import Counter, defaultdict
import os

class PositionState:
	nucl = ["A", "G", "C", "T", "-"]
	# Match emission probability
		# {"A":0.5, "G": 0.7, "C":0.1, "T":0.6}

	# Delete emission probability
	delete_dict = {"A":0, "G": 0, "C":0, "T":0, "-":1}

	# Insertion emission probability, background prob of nucleotide
		# {"A:0.2"}
	def __init__(self, bases,PFM, pos):
		self.pos = pos
		match_dict = create_match(bases,PFM)
		insert_dict = create_insert(bases,PFM)

	def create_match(bases,PFM):
		match_dict = {}
		num_seqs = len(bases)
		num_gaps = bases.count("_")

		if num_gaps >= num_seqs/2: 	# match state includes only colums where there are bases in at least half the sequences (no gaps)
			match_dict = None
		else:
			for idx in range(4):
				match_dict[nucl[idx]] = PFM[nucl[idx]][pos]/num_seqs # "add-one rule"

		return match_dict

	# insert is the background probability of nucleotide (in whole genome or the conserved region??)	
	def create_insert(bases,PFM):
		insert_dict = {}
		length_seq = PFM.shape()[1]
		for idx in range(4):
			insert_dict[nucl[idx]] = PFM[idx].sum()/length_seq

		return insert_dict



class pHMM:
	# transition probabilities btw states
		# 3x3 Matrix
	# emission probabilities by position

	def __init__(self, sequences):
		self.transition = create_transition(sequences)
		self.emission = create_emission(sequences)

	def create_transition(sequences):
		pass

	def create_emission(sequences):
		emission = []
		for pos,bases in enumerate(zip(*sequences)): # bases = tuple with all the letters at the position in it
			emission.append(PositionState(bases,PFM,pos))

		return emission

	def get_PFM(sequences):
		PFM = np.zeros((4,len(sequences[0])))
		for pos,bases in enumerate(zip(*sequences)): # list of tuples w/ 0th tuple being all the letters in 0th position, etc
			PFM[0][pos] = bases.count("A")
			PFM[1][pos] = bases.count("G")
			PFM[2][pos] = bases.count("C")
			PFM[3][pos] = bases.count("T")
			PFM[4][pos] = bases.count("-")

		return PFM


class GeneSequence:
	def __init__(self, iterator, genome):
		self.sequence = ""
		self.gene_name = None
		self.location = None
		self.genome = genome

		for match in iterator:
			text = match.group()
			print(text)
			if "gene=" in text:
				self.gene_name = text[5:]
			if "location=" in text:
				 #get rid of non-numeric and non-periods
				self.location = re.sub(r"[^0-9\.]",r"", text).split("..")

def read_in_annotated(filename, p): # p is the compiled pattern object
	current_gene_sequence = None
	gene_sequences = []

	with open(filename,'r') as f:
		print("opened file")
		count = 0
		for line in f:
			if ">lcl" in line and "gene=" in line: # start of new annotated
				iterator = p.finditer(line)
				genome = filename.split("/")[-1]
				current_gene_sequence = GeneSequence(iterator,genome.split("_")[0])
				gene_sequences.append(current_gene_sequence)
			elif ">lcl" in line and "gene=" not in line:
				current_gene_sequence = None

			if current_gene_sequence is not None and ">lcl" not in line:
				current_gene_sequence.sequence += line[:-1] #don't include \n
				

	# Return a list of the genes surrounding dmDa for a genome
	return gene_sequences

def create_gene_dict(gene_list):
	gene_dict = defaultdict(list)
	for gene in gene_list:
		gene_dict[gene.gene_name].append(gene)
	# input = list of all sequences
	# output = a dict where the key is the name of gene and the val is a list of all the gene sequences

	#maybe remove the genes that only appear once
	return gene_dict

def create_pHMMs(gene_dict):
	# input = dict from above function
	# for each key-val pair, align the sequences and create a pHMM
	# output = dict where name of gene is the key, and the val is the pHMM
	pHMM_dict = defaultdict(list)
	for k,v in gene_dict.items():
		pHMM_dict[k].append(pHMM(v))

	return pHMM_dict

def score_sequence(region_seq): # viterbi or forward
	# input = conserved region
	# get a score using every pHMM
	# output = a dict with gene name as key, and val is the score
	pass


def create_unsorted_gene_list():
	path = "data/annotations/"
	filenames = [file for file in next(os.walk(path))[2] if ".txt" in file] #["Dshibae_annotated.txt", "Oantarcticus307_annotated.txt", "Oarcticus238_annotated.txt", "Rmobilis_annotated.txt", "Rpomeroyi_annotated.txt"]
	
	unsorted_gene_list = []

	p = re.compile("(?<=\[).+?(?=\])") #regex expression to match things in brackets
	# read all files in annotation folder
	for file in filenames:
		gene_list = read_in_annotated(path+file,p)
		unsorted_gene_list.extend(gene_list)

	# get total gene list to worry about
	# output = unsorted gene list
	return unsorted_gene_list

'''
unsorted_gene_list = create_unsorted_gene_list()
gene_dict = create_gene_dict(unsorted_gene_list)
'''

# Create profile HMM model of a gene
	# list - gene = [seq1, seq2, seq3]
	# get their alignments
	# get modified PFM for all the sequences (to check whether a col is mostly gaps)
	# states = insertion, deletion, match
	# match state (A,G,C,T) = freq of base at that position/# seq, add a dirichlet for the 0s
	# transition probabilities = number of transitions from i to j/number of transitinos from i to other states

# Create profile HMM model for all genes

# Get probability of gene name based on sequence using profile HMM

# Get gene name from BLAST and get e-value

# Combine and predict a gene name