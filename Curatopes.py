#!/usr/bin/env python

import numpy as np
import pandas as pd
import math
import sys, glob, os, re
import getopt
from functools import reduce
import subprocess
from lib import CuraLib
scriptPath = os.path.abspath(os.path.dirname(__file__))

# What this will do:
# Read FPKM from samples
# Read GTEx statistics
# Return: predictions of general population-applicable epitopes
# Test Command:
# ./Curatopes.py -i "./Data/Samples/TestSet.115821.NoV"  -o "./Data/" -e "Melanoma"

helpMessage="""\
Curatopes.py
	-e <Tumor entity>\te.g Melanoma, Uveal Melanoma, Glioblastoma, etc.
	-i <inputfile>\tExpression tumorData generated from biosamples with gene x sample format [TPM / FPKM] - Headerless - ENSG-Versionless
	-o <output dir>\tWill make a directory for a given entity to store output
	-m <MHC alleles>\tProvide a list of MHC
"""
#<<< <<< <<< Hardcoded files >>> >>> >>>#
# Use HPA low+medium+high
hpaList = "../Data/HPA/high+medium+low"
# GTEx tissue files need to be parsed separately
TissueStats = '../Data/GTEX/TissueStats'
# Proteome -
# ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
ProteomeAnnotation = "../Data/Proteome/Homo_sapiens.GRCh38.pep.all.protein_coding.txt"
# Critical Tissue List
CriticalTissue = "../Data/CriticalTissue/criticaltissues_binary.csv"
# list of genes from version 1 for comparison
cmV1 = "../Data/curatopes-melanoma-v1-genes.txt"
# network interaction scores
NIS = "../Data/Curatopes_UvealMelanoma/network-interaction-scores.csv"

#<<< <<< <<< Hardcoded End >>> >>> >>>#

def status(*message):
	print(*message, file = sys.stderr)


def removeVersion(x):

	return x.split(".")[0]


def main(argv):

	inputfile = ''
	outputfile = ''
	mhcfile = ''
	entity = ''

	try:

		opts, args = getopt.getopt(argv, "hi:m:o:e:", ["inputfile=", "outputfile=", "entity=", "mhcfile="])

	except getopt.GetoptError:

		status(helpMessage)
		sys.exit(2)

	for opt, arg in opts:

		if opt == '-h':

			status(helpMessage)
			sys.exit()

		elif opt in ("-i", "--input"):

			inputfile = arg
			status("Expression Input is", inputfile)

		elif opt in ("-m", "--mhc"):
			mhcfile = arg
			status("MHC profile file is", mhcfile)

		elif opt in ("-o", "--output"):

			outputfile = arg
			status("Output file is", outputfile)

		elif opt in ("-e", "--entity"):

			entity= arg
			status("Tumor entity is", entity)


	return [inputfile, outputfile, entity, mhcfile]

if __name__ == "__main__":
	inFile, outFile, tumorEntity, mhcAlleles = main(sys.argv[1:])

if len(sys.argv[1:]) < 1:
	sys.exit(status(helpMessage))


# read inputfile, no header, genenames on the rows
tumorData = pd.read_table(inFile, header = None, index_col = 0)

# aggregate redundant ENSGs (e.g. ENSG00000002586_PAR_Y and ENSG00000002586) into one row (by summing expression values)
tumorData = tumorData.groupby(tumorData.index.map(removeVersion), axis = 0, as_index = True).sum()

hpaFilter = pd.read_table(hpaList, header = None, index_col = 0)
hpaFilter.index = hpaFilter.index.map(removeVersion)

with open(mhcAlleles, "r") as f:
    hla = [line.strip() for line in f]

status("(Re-)applying TPM conversion...")
tumorData = tumorData / tumorData.sum(axis = 0) * 1e6
print(tumorData.head())

# Remove genes with row sums < ten
status("Processing expression data...")

# Calculate summary statistics for the tumor.
### tumorData['TumorQ10'] = tumorData.quantile(.1, axis = 1)

### HIEGE MOD -> Increase tumor quantile to 30th
tumorData['TumorQ10'] = tumorData.quantile(.30, axis = 1)

tumorData['TumorMedian'] = tumorData.median(axis = 1)

filter_labels = ["TuQ10gt1", "notInHPA", "ProtAnnot", "inTissue", "TuQ10gtTiQ90", "CritTissues"]
filters = pd.DataFrame({f"f{d}_{n}": [None] * len(tumorData.index) for d, n in enumerate(filter_labels)})
filters.index = tumorData.index

filters.iloc[:, 0] = tumorData['TumorQ10'] > 1

# intersect HPAFilter with tumorData
status("HPA filtering...")
filters.iloc[:, 1] = ~tumorData.index.isin(hpaFilter.index)

# filter for protein-coding genes only
status("Proteome annotation...")
proteomeData = pd.read_table(ProteomeAnnotation, header = None, index_col = 1)
filters.iloc[:, 2] = tumorData.index.isin(proteomeData.index)

# Preparing for GTEx data parsing
all_tissues = glob.glob(TissueStats + "/*.stat")

tissue_list = []

status("Tissue annotation...")
for filename in all_tissues:

	tissueName = os.path.basename(filename).split(".")[0]
	tissues_stats = pd.read_table(filename, index_col = 0, header = 0)
	tissues_stats.index = tissues_stats.index.map(removeVersion)

	tissues_stats.rename(columns = {"qu90": tissueName}, inplace = True)

	tissue_list.append(tissues_stats[[tissueName]])

status(f"Processed {len(tissue_list)} tissues")


# Merging all tissue q90s
data_tpm = reduce(lambda x, y: pd.merge(x, y, left_index = True, right_index = True, how = "outer"), tissue_list)
filters.iloc[:, 3] = tumorData.index.isin(data_tpm.index)


# get genes for which we have both tumor and tissue data
h = tumorData.index.intersection(data_tpm.index)

# filter genes whose tumor expression q10 is consistently higher than their tissue expression q90
data_bin = data_tpm.loc[h, :].lt(tumorData.loc[h, 'TumorQ10'], axis = 0)
filters.loc[h, filters.columns[4]] = data_bin.sum(axis = 1) / data_bin.shape[1]

# load the GTEx subset of critical tissues
critList = pd.read_table(CriticalTissue, header = 0, index_col = 0)
critList = critList[critList.selected == 1].index

# extract expression of all critical tissues
data_crit = data_tpm.loc[:, data_tpm.columns.isin(critList)]
data_crit = data_crit < 10
filters.loc[h, filters.columns[5]] = data_crit.loc[h, :].sum(axis = 1) / data_crit.shape[1]

# output filter result table
outPath = os.path.join(outFile, tumorEntity)
os.makedirs(outPath, exist_ok = True)

filters.to_csv(os.path.join(outPath, "filters.csv"), sep = "\t")

### BENCHMARKING BLOCK ###
# # output expression table for genes selected in version 1 of Curatopes Melanoma
# with open(cmV1) as f:
#     V1 = [l.strip() for l in f if l[0] != "#"]
#     tumorData.loc[V1, ["TumorQ10"]].merge(
#         data_tpm.loc[V1, :], left_index = True, right_index = True
#     ).to_csv(os.path.join(outPath, "V1.tpms"), sep = "\t")


# identify and output tier 1 genes
tier1 = filters.index[(filters == 1).apply(sum, axis = 1) == filters.shape[1]]
with open(os.path.join(outPath, "tier1.txt"), "w") as f:
	f.write("\n".join(tier1.tolist()))
# identify and output tier 2 genes
tier2 = filters.index[(filters.iloc[:, :-1] == 1).apply(sum, axis = 1) == filters.shape[1] - 1]
tier2 = tier2.difference(tier1)
with open(os.path.join(outPath, "tier2.txt"), "w") as f:
	f.write("\n".join(tier2.tolist()))
### BENCHMARKING BLOCK ###



# run predictions of epitopes and their parameters; this is the time-consuming step
status("Performing predictions...")
proteomeDict = CuraLib.createProteomeDict()
data_pred = CuraLib.runPreds(HLAList = hla, GeneList = (tier1.tolist() + tier2.tolist()), Proteome = proteomeDict, ic50cutoff = None)
status("Predictions complete...")

# convert result to DataFrame
data_pred = pd.DataFrame.from_records(data_pred)
data_pred = data_pred.set_index(data_pred[0])

# calculate the fourth coefficient (AKA "expression index") for the gPIE index
status("Calculating f4...")
upperBound = data_tpm.max(axis = 1)
data_tpm.loc[h, 'f4'] = tumorData.loc[h, 'TumorQ10'] / (tumorData.loc[h, 'TumorQ10'] + upperBound[h])  # TODO: keep f4 separate

# Merge expression data with predicton data
status("Merging data for output...")
data_comb = data_pred.merge(data_tpm['f4'], left_index = True, right_index = True, how = "left")
data_comb = data_comb.merge(tumorData['TumorMedian'], left_index = True, right_index = True, how = "left")

status("Computing physiochemical properties...")
Physchem = pd.Series(data_comb[3].unique()).apply(CuraLib.calcImmunochem)

# Headers
Physchem = pd.DataFrame(Physchem.to_list(),
	columns = [
	"Peptide",
	"Iso",
	"MolW",
	"Stability",
	"Hydrophobicity",
	"Polarity"])

data_comb.columns = ['GeneID', 'TranscriptID', 'HLA-Allele', 'Peptide', 'IC50', 'Immunogenicity', 'ExpressionIndex', 'TumorMedian']

data_comb = data_comb.merge(Physchem, on='Peptide')



status("Computing new gPIE...")
#	1	2	3	4	5	f4	TumorMedian	Iso	MolW	Stability	Hydrophobicity	Polarity
#0	ENST00000443026	HLA-A0101	VTDEEMERA	97.35	0.15489	0.8246952519698	53.04549019606092	4.0500284194946286	1079.1392	35.67777777777778	-1.2555555555555555	28.29111111111111
#data_comb.columns = ['TranscriptID','HLA-Allele','Peptide','IC50','Immunogenicity','ExpressionIndex','TumorMedian','Iso','MolW','Stability','Hydrophobicity','Polarity']
epitopeFile = os.path.join(outPath, "Results.csv.gz")
data_comb.to_csv(epitopeFile, sep = "\t", compression = "gzip")


# Retrieve path to Rscript
PathToR = subprocess.run(['/bin/bash', '-c', 'which Rscript'],
	stdout = subprocess.PIPE,
	universal_newlines = True)

PathToR = PathToR.stdout.rstrip()

# Compute activity and binding model
try:
	Score = subprocess.call([PathToR, os.path.join(scriptPath, 'lib', 'RF-ModelsCuratopes.R'), epitopeFile, NIS, outPath])
except Exception as e:
	sys.exit("Scoring not working.")

status("Everything is fine! Curatopes pipeline is done.")
