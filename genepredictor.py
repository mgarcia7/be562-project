################################################################
# TASK: Predict what gene it is based on annotated genes + BLAST
# INPUT: Conserved regions
# OUTPUT: Name of gene of all the conserved regions/likelihood that it is that gene
################################################################
import numpy as np
import re

class PositionState:
	# Match emission probability
	# Delete emission probability
	# Insertion emission probability

	def __init__(self, sequences):
		self.match = create_match(sequences)
		self.delete = create_delete(sequences)
		self.insert = create_insert(sequences)

	def create_match(sequences):
		pass

	def create_delete(sequences):
		pass

	def create_insert(sequences):
		pass


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
		pass


class GeneSequence:
	def __init__(self, iterator, genome):
		self.sequence = ""
		self.gene_name = None
		self.location = None
		self.genome = genome

		for match in iterator:
			text = match.group()
			if "gene=" in text:
				gene_name = text[5:]
			if "location=" in text:
				location = re.sub('[^0-9]\.','', a) #get rid of non-numeric and non-periods
				location = tuple(text.split(".."))


def read_in_annotated(filename, p): # p is the compiled pattern object
	# Read in line by line, save sequence, check every 5 if the 5th one is dmDA
		# if not, clear list and start over
		# if so, save list and do the next 5 lines after dmDa
	current_gene_sequence = None
	gene_sequences = []

	with open(filename,'r') as f:
		count = 0
		for line in f:
			if ">lcl" in line and "gene=" in line: # start of new annotated
				iterator = p.finditer(line)
				current_gene_sequence = GeneSequence(iterator,filename.split("_")[0])

			else if ">lcl" in line and "gene=" not in line:
				current_gene_sequence = None

			if current_gene_sequence is not None:
				current_gene_sequence.sequence += line
				



	# Return a list of the genes surrounding dmDa for a genome

	pass

def create_gene_dict(genome_sequences)
	# input = list of all sequences
	# sort them by gene name and add them to a dict
	# output = a dict where the key is the name of gene and the val is a list of all the gene sequences
	pass

def create_pHMMs(gene_sequences):
	# input = dict from above function
	# for each key-val pair, align the sequences and create a pHMM
	# output = dict where name of gene is the key, and the val is the pHMM
	pass

def score_sequence(region_seq):
	# input = conserved region
	# get a score using every pHMM
	# output = a dict with gene name as key, and val is the score
	pass


def create_unsorted_gene_list():
	path = "data/annotations/"
	filenames = ["Dshibae_annotated.txt", "Oantarcticus307_annotated.txt", "Oarcticus238_annotated.txt", "Rmobilis_annotated.txt", "Rpomeroyi_annotated.txt"]
	unsorted_gene_list = []

	p = re.compile("(?<=\[).+?(?=\])") #regex expression to match things in brackets
	# read all files in annotation folder
	for file in filenames:
		gene_list = read_in_annotated(path+file,p)
		unsorted_gene_list.append(gene_list)

	# get total gene list to worry about
	# output = unsorted gene list
	return unsorted_gene_list

# Create profile HMM model of a gene
	# list - gene = [seq1, seq2, seq3]
	# get their alignments
	# get modified PFM for all the sequences (to check whether a col is mostly gaps)
	# states = insertion, deletion, match
	# match state (A,G,C,T) = freq of base at that position/# seq, add a dirichlet for the 0s
	# transition probabilities = number of transitions from i to j/number of transitinos from i to other states

# Create profile HMM model for all genes

# Get probability of gene name based on sequence using profile HMM

# Get gene name from BLAST and get e value

# Combine and predict a gene name