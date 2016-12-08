# This program creates a csv file with the names of the organism, and starting and ending position

import genepredictor as gp
from collections import Counter, defaultdict

unsorted_gene_list = gp.create_unsorted_gene_list()
#gene_dict = gp.create_gene_dict(unsorted_gene_list)

filenames = ["Dshibae", "Oantarcticus307", "Oarcticus238", "Rmobilis", "Rpomeroyi", "Roseovarius217", "Rlitoralis", "Rdentrificans", "JannaschiaCCS1", "Lvestfoldensis", "Rbacterium"]

outfile = open('transformed.csv', 'w')

gene_names = ["dmda", "dddp","dddl","dddd","dddq","dddw"]
rows = defaultdict(list)
for file in filenames:
	rows[file] = ["x", "x", "x", "x", "x", "x"]

# Name of Genome	dmdA	dddP	dddL	dddD	dddQ	dddW
for k in unsorted_gene_list:
	if k.gene_name.lower() in gene_names:
		ind = gene_names.index(k.gene_name.lower())
		rows[k.genome][ind] = "..".join(k.location)


# print headings
header_line = ','.join(["genome","dmdA", "dddP","dddL","dddD","dddQ","dddW"])
print(header_line, file=outfile)

# goes through dict and prints info to file
for key,val in rows.items():
	info = [key]
	info.extend(val)
	newline = ",".join(info)
	print(newline, file=outfile)

		

outfile.close()



