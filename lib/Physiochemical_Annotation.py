#!/usr/bin/python3

import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis

### Files needed are the List of Epitopes.

pep = sys.stdin.read().split()

input_allele = sys.argv[1]

for input_pep in pep:

	analysed_seq = ProteinAnalysis(input_pep)

	molW = analysed_seq.molecular_weight()

	stability = analysed_seq.instability_index()

	iso = analysed_seq.isoelectric_point()

	hyrdophobicity = analysed_seq.gravy()

	polarity = {'A':0.00, 'R':52.00, 'D':49.70, 'N':3.38, 'C': 1.48,'E':49.90, 'Q':3.53, 'G':0.00, 'H':51.60, 'L':0.13, 'I' : 0.13, 'K':49.50, 'M':1.43, 'F':0.35, 'P':1.58, 'S':1.67, 'T':1.66, 'W':2.10,'Y':1.61, 'V':0.13}

	sumPol = 0

	for i in input_pep:
		
		if i in polarity:

			sumPol += polarity[i]

	sumPol = sumPol/len(input_pep)

	##### allele #### peptide #### IP #### molW #### stability #### gravy #### polarity  #### ID 
	print(input_allele,'\t',input_pep,'\t',iso,'\t',molW,'\t',stability,'\t',hyrdophobicity,'\t',sumPol)

### POLARITY VALUES HARDCODE OF DOOM ###
# SOURCE: 
# J. Theoret. Biol. (1968) 21, 170-201
# The Characterization of Amino Acid Sequences in Proteins by Statistical Methods 
# Zimmerman, J.M et al.
# DOI: 10.1016/0022-5193(68)90069-6

# Ala 	A	0.00
# Arg 	R	52.00
# Asp	D	49.70
# Asn 	N	3.38
# Cys 	C	1.48
# Glu	E	49.90
# Gln	Q	3.53
# Gly	G	0.00
# His	H	51.60
# Leu	L	0.13
# Ile 	I	0.13
# Lys 	K	49.50
# Met	M	1.43
# Phe	F	0.35
# Pro 	P	1.58
# Ser	S	1.67
# Thr 	T	1.66
# Trp	W	2.10
# Tyr 	Y	1.61
# Val 	V	0.13