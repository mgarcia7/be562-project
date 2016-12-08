################################################################
# TASK: Predict what gene it is based on annotated genes + BLAST
# INPUT: Conserved regions
# OUTPUT: Name of gene of all the conserved regions/likelihood that it is that gene
################################################################
import numpy as np
import re
from collections import Counter, defaultdict
import os
from math import log, floor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import Bio.SeqIO as SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
from itertools import cycle
import time

class PositionState:
	nucl = ["A", "G", "C", "T", "-"]
	# Match emission probability
		# {"A":0.5, "G": 0.7, "C":0.1, "T":0.6}

	# Delete emission probability
	delete_dict = {"A":0, "G": 0, "C":0, "T":0, "-":1}

	# Insertion emission probability, background prob of nucleotide
		# {"A:0.2"}
	def __init__(self, bases,PFM, pos,len_seq):
		self.len_seq = len_seq
		self.pos = pos
		self.match_dict = self.create_match(bases,PFM)
		self.insert_dict = self.create_insert(bases,PFM)

	def create_match(self,bases,PFM):
		match_dict = {}
		num_seqs = len(bases)
		for idx in range(4):
			match_dict[self.nucl[idx]] = (PFM[idx][self.pos]+1)/(num_seqs+4) # implement the "add-one rule"

		return match_dict

	# insert is the background probability of nucleotide (in whole genome or the conserved region??)	
	def create_insert(self,bases,PFM):
		insert_dict = {}
		num_seqs = len(bases)


		for idx in range(4):
			insert_dict[self.nucl[idx]] = (PFM[idx][self.pos].sum()+1)/(num_seqs + 4)

		return insert_dict

class pHMM:
	# match state as defined
	# insert state if not match state
	# delete state if it's technically a match or insert state but that particular sequence has a gap
	# transition probabilities btw states
		# 3x3 Matrix
			# M D I
		# M
		# D
		# I
	# emission probabilities by position

	def __init__(self, sequences):
		self.PFM = get_PFM(sequences)
		self.emission, self.states = self.create_emission(sequences)
		self.transition = self.create_transition(sequences)

	def create_transition(self, sequences):
		state_dict = {"M":0, "D":1, "I":2}
		transition = np.zeros((3,3))
	
		for seq in sequences:
			current_states_list = self.states[:] #copy of states
			for pos,base in enumerate(seq): #assume state -1 is a match
				if pos == 0:
					prev = "M"
				else:
					if base == "-":
						current_states_list[pos] = "D"

					prev = current_states_list[pos-1]
					current = current_states_list[pos]

					if prev == "I" and prev == current and (pos == 0 or pos == None):
						transition[state_dict[prev]][state_dict[current]] += 1.5
					else:
						transition[state_dict[prev]][state_dict[current]] += 1


		return (transition+1)/(transition.sum(axis=0)+3) # returns normalized transition matrix by col
 
	def create_emission(self, sequences):
		emission = []
		states = []
		for pos,bases in enumerate(zip(*sequences)): # bases = tuple with all the letters at the position in it
			num_gaps = 0
			for base in bases:
				if base == "-":
					num_gaps+=1

			if num_gaps >= floor(len(sequences)/2.0): #beginning state, end state, match states need to not have many gaps
				states.append("I")
				continue
			else:
				states.append("M")
				emission.append(PositionState(bases,self.PFM,pos,len(sequences[0])))

		#emission.append(None) # signals end state
		return emission,states

def get_PFM(sequences):
	PFM = np.zeros((5,len(sequences[0])))
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
			if "gene=" in text:
				self.gene_name = text[5:]
			if "location=" in text:
				self.location = re.sub(r"[^0-9\.]",r"", text).split("..") #get rid of non-numeric and non-periods

def read_in_annotated(filename, pattern): # p is the compiled pattern object
	current_gene_sequence = None
	gene_sequences = []

	with open(filename,'r') as f:
		count = 0
		for line in f:
			if ">lcl" in line and "gene=" in line: # start of new annotated
				iterator = pattern.finditer(line)
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
	pHMM_dict = dict()

	for k,v in gene_dict.items():
		if len(v) >= 14:
			aligned = multiseqalign(k,v)
			pHMM_dict[k] = pHMM(aligned)
	
	#aligned = ["AG---C","AGAG-C", "AG-AAG", "-GAAAC", "AG---C"]
	#pHMM_dict["test"] = pHMM(aligned)

	return pHMM_dict

def score_sequence(region_seq,pHMM_dict): # viterbi
	# input = conserved region
	# get a score using every pHMM
	score_dict = dict()
	for gene_name, pHMM in pHMM_dict.items():
		print("Scoring w", gene_name, "model...")
		score_dict[gene_name] = viterbi(pHMM,region_seq)
		if score_dict[gene_name] == None:
			del score_dict[gene_name]

		#print("Done scoring")
	# output = a dict with gene name as key, and val is the score
	return score_dict


def viterbi(pHMM,sequence):
	L = pHMM.states.count("M") # number of matching states
	N = len(sequence) #length of seq
	em = cycle(pHMM.emission)
	tr = pHMM.transition

	#try:
	Fm = np.zeros((N,L))
	Fi = np.zeros((N,L))
	Fd = np.zeros((N,L))

	Fm[0][0] = 1 #initialize

	current_state_em = next(em)

	for j in range(N): # goes through sequence
		if j == 0:
			continue
		base = sequence[j]
		for i in range(L):
			Fm[j][i] = log(current_state_em.match_dict[base]) + max(Fm[j-1][i-1]+log(tr[0][0]), Fi[j-1][i-1] + log(tr[2][0]), Fi[j-1][i-1] + log(tr[1][0]))
			Fi[j][i] = log(current_state_em.insert_dict[base]) + max(Fm[j][i-1]+log(tr[0][2]), Fi[j][i-1]+log(tr[2][2]))
			Fd[j][i] = max(Fm[j-1][i]+log(tr[0][1]),Fd[j-1][i]+log(tr[1][1]))
			current_state_em = next(em)

			if max(Fm[j][i],Fi[j][i],Fd[j][i]) < -5:
				return None

	return max(Fm[N-1][L-1],Fi[N-1][L-1],Fd[N-1][L-1])
	#except:
	#	return None

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

def multiseqalign(gene_name,sequences):
	in_file = "data/alignments/unaligned.fasta"
	out_file = "data/alignments/aligned" + "_" + gene_name + ".fasta"

	if not os.path.isfile(out_file):
		print("Aligning ", gene_name, "....")
		# Write my sequences to a fasta file
		records = (SeqRecord(Seq(gene.sequence,generic_dna), id=str(index), name="Test", description="Test") for index,gene in enumerate(sequences) )

		with open(in_file, 'w') as handle:
			SeqIO.write(records, handle, "fasta")

		clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, auto=True, force=True)
		clustalomega_cline()
		print("Aligned")

	alignment = AlignIO.read(out_file,format="fasta")
	return [str(record.seq) for record in alignment]

def initialize_pHMM_models():
	unsorted_gene_list = create_unsorted_gene_list()
	gene_dict = create_gene_dict(unsorted_gene_list)
	phmm_dict = create_pHMMs(gene_dict)
	return phmm_dict

def score_all_conserved_regions(conserved_regions,pHMM_dicts):
	all_score_dicts = []
	for region in conserved_regions:
		c = region.conserved.region1.seq
		c = c.replace("-","")
		all_score_dicts.append(score_sequence(c))

	# returns list of dicts
	return all_score_dicts

def predictions(conserved_regions,all_score_dicts):
	for idx,score_dict in enumerate(all_score_dicts):
		for i in range(2):
			maxkey = keywithmaxval(score_dict)
			maxscore = score_dict[keywithmaxval(score_dict)]
			conserved_regions[idx].possible_genes.append((maxkey,maxscore))
			del score_dict[maxkey]
	return conserved_regions

def print(conserved_regions):
	for regions in conserved_regions:
		print(regions.conserved.region1.organism,regions.conserved.region2.organism)
		print(regions.conserved.possible_genes)
		print()



def keywithmaxval(d):
     """ a) create a list of the dict's keys and values; 
         b) return the key with the max value"""  
     v=list(d.values())
     k=list(d.keys())
     return k[v.index(max(v))]

'''

start_time = time.time()

unsorted_gene_list = create_unsorted_gene_list()
gene_dict = create_gene_dict(unsorted_gene_list)
phmm_dict = create_pHMMs(gene_dict)

possible_gene = []

for i in range(2):
	maxkey = keywithmaxval(score_dict)
	maxscore = score_dict[keywithmaxval(score_dict)]
	possible_gene.append((maxkey,maxscore))
	del score_dict[maxkey]

print(possible_gene)
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