################################################################
# TASK: Predict what gene it is based on annotated genes + BLAST
# INPUT: Conserved regions
# OUTPUT: Name of gene of all the conserved regions/likelihood that it is that gene
################################################################
import numpy as np

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
		self.sequences = sequences
		self.transition = create_transition(sequences)
		self.emission = create_emission(sequences)

	def create_transition(sequences):
		pass

	def create_emission(sequences):
		pass


def 

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