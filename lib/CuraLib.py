#!/usr/bin/env python

import tempfile
import sys, glob, os, re
import subprocess

from lib import predict_immunogenicity_mod as predImmu
import concurrent.futures
import multiprocessing as mp
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# required tools from commons - take care of that later when packaging it.
PRED="/home/_common/Apps/IEDB/mhc_i/src/predict_binding.py"
NETMHC="/home/_common/Apps/NetMHCpan/netMHC-4.0/netMHC"

# Return all natively providd HLA alleles in netMHCpan
# hlaList can hold all HLA alleles that we can predict from with default settings
NETMHC_HLAS = subprocess.run([NETMHC, '-listMHC'], stdout = subprocess.PIPE, universal_newlines = True)
NETMHC_HLAS = set(re.findall(r"^HLA-\S+", NETMHC_HLAS.stdout, re.MULTILINE))

# Accept lists of genes and aquire proteome sequences.
# Return Peptide - IC50 - GeneID
# Test list of Genes
# Hardcode - Only change when updated
Proteome = "../Data/Proteome/Homo_sapiens.GRCh38.pep.all.linear.fa"

# Perl code to linerarize a fasta file - just do it.
# perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' in.fasta > out.fasta
# Input a proteome fasta file linearized
# This provides a path to a folder containing all relevant custom hla alleles
hlaFolder = "../Data/HLA_Sequences"

CUSTOM_HLAS = dict()
for fname in glob.glob(os.path.join(hlaFolder, "*-protein.fa")):
        label = os.path.basename(fname).replace("-protein.fa", "")
        CUSTOM_HLAS[label] = fname


def createProteomeDict(Proteome = Proteome):

        proteomeDict = defaultdict(lambda: defaultdict(dict))
        ensembl_pattern = re.compile(r"ENS.[0-9]{11}")

        for record in SeqIO.parse(Proteome, "fasta"):

                # sorted order is ENSG..., ENSP..., ENST...
                gene, protein, transcript = sorted(ensembl_pattern.findall(record.description))

                proteomeDict[gene][transcript] = (protein, str(record.seq))

        return proteomeDict


# Create a subset of the proteome either with interseciton option
# complement = F: Create a proteome list of all the genes of interest
# complement = T: Remove all Protein sequences belonging to a particular list of genes of interest from the proteome
# Flatten in this case returns a list of values
def createProteomeSlice(proteomeDict, GeneList, complement = False, flatten = False, pepslice = False):

        if flatten and pepslice:

                raise ValueError("Cannot both flatten and pepslice!")

        if not complement:

                proteomeSlice = {gene : proteomeDict[gene] for gene in GeneList & proteomeDict.keys()}

        else:

                proteomeSlice = {gene : proteomeDict[gene] for gene in proteomeDict.keys() - GeneList}


        if flatten:

                flatList = []

                for gene in proteomeSlice:
                        for transcript in proteomeSlice[gene]:
                                flatList.append(proteomeSlice[gene][transcript][1])

                return "\n".join(flatList)

        elif pepslice:

                #d = defaultdict(int)
                d = set()

                for geneseqs in proteomeSlice.values():
                        for peptide, seq in geneseqs.values():
                                for l in range(9, 13):
                                        for p in range(len(seq) - l + 1):
                                                #d[seq[p : p+l]] += 1
                                                d.add(seq[p : p+l])

                return d

        else:

                return proteomeSlice


# Filter the epitome against the proteome
def redundantSequence(Peptide, Proteome):

        return Peptide in Proteome


###################################################################
###################################################################
# Predictor functions should return identifiable lists of predictions
# Functions should be desgined in such a way that we can easily add new features.
def predictIC50(ID, Allele, ProteinSequence, Proteome = ''):

        # NOTE: NetMHCpan does not accept sequence input from stdin, so we have to work around that with temporary files.

        # prepare protein sequence in FASTA format
        fileContent = f">{ID}\n{ProteinSequence}"

        tempFasta = tempfile.NamedTemporaryFile(delete = True, suffix = '.fa')
        tempFasta.write(fileContent.encode('ascii'))
        tempFasta.flush()


        pepList = []

        # run netMHCpan predictions
        # TODO: We need to compensate for unavailable MHC alleles.

        if Allele in NETMHC_HLAS:

                prediction = subprocess.run([NETMHC, '-s', '-l', '9,10,11,12', '-a', Allele, '-f', tempFasta.name],
                        stdout = subprocess.PIPE, universal_newlines = True)
                pred = prediction.stdout

                # This excludes the netMHCpan header.
                predList = pred.split("\n")[42:]

                for elements in predList:

                        # Processing of netMHCpan's output
                        # This excludes the netMHCpan footer.
                        if not elements.startswith("Protein"):

                                ele = re.split(r"\s+", elements)

                                try:
                                        # Remove gene ID and HLA allele from result since we can identify it through the proteome
                                        line = [ele[i] for i in (11, 2, 3, 13)]

                                except IndexError as e:

                                        continue

                                else:

                                        # Remove "X" peptides
                                        if not "X" in line[2]:

                                                # ['ENST00000335137', 'HLA-A0201', 'KDMKTAIRQLRK', '35012.54']
                                                if not redundantSequence(line[2], Proteome):
                                                            pepList.append(line)

        elif Allele in CUSTOM_HLAS:

                # predict_binding_modified.py netmhcpan -m $i $j "$PROTS" > "${TOPE_FOLDER}/${NAME}/EpiDB.${NAME}.${j}l.raw"

                for predLength in range(9, 13):

                        prediction = subprocess.run(['python2', PRED, 'netmhcpan', '-m', CUSTOM_HLAS[Allele], str(predLength), tempFasta.name],
                                stdout = subprocess.PIPE, universal_newlines = True)
                        pred = prediction.stdout

                        predList = pred.split("\n")

                        for elements in predList:

                                if not elements.startswith("allele"):

                                        ele = elements.split("\t")

                                        if len(ele) == 4:

                                                ele[0] = ID

                                                if not "X" in ele[2]:

                                                        if not redundantSequence(ele[2], Proteome):

                                                                ele[1] = Allele

                                                                pepList.append(ele)

        else:

                raise FileNotFoundError(f"Custom allele {Allele} not found in {hlaFolder} for use with netMHCpan!")

        tempFasta.close()  # make sure that the file is deleted

        return pepList


def predictImmu(Peptide, Allele, ID = ''):

        prediction = predImmu.Prediction()

        x = prediction.predict((Peptide, None, Allele))

        x = list(x)

        return x


# def predictTap
# def predictPhyschem
###################################################################
###################################################################
#Input: HLA-List,Fasta of genes,or hash slice of proteome.
def runPreds(HLAList, GeneList, Proteome, ic50cutoff = 500):

        if type(ic50cutoff) not in (int, float) and ic50cutoff != None:
                raise TypeError(f"Invalid value for parameter ic50cutoff: {ic50cutoff}")

        # Create prediction hla allele list
        predSet = HLAList
        # Create a hash slice for proteins in question
        proteomeSlice = createProteomeSlice(Proteome, GeneList)
        print(f"Created protein sequences for {len(proteomeSlice)} genes", file = sys.stderr)

        # Create set difference proteome
        #setDiffProteome = createProteomeSlice(ProteomeDict, GeneList, complement = True, pepslice = True)
        setDiffProteome = createProteomeSlice(Proteome, GeneList, complement = True, flatten = True)
        print(f"Created peptide slice of length {len(setDiffProteome)}", file = sys.stderr)

        # Storage for the returned epitopes
        epitome = []
        # Processes list
        transcriptsToGeneMap = dict()

        with concurrent.futures.ProcessPoolExecutor(max_workers = 60) as executor:

                pred = []

                for allele in predSet:
                        for gene in proteomeSlice:
                                for transcript in proteomeSlice[gene]:

                                        transcriptsToGeneMap[transcript] = gene

                                        result = executor.submit(predictIC50, ID = transcript, Allele = allele,
                                                ProteinSequence = proteomeSlice[gene][transcript][1], Proteome = setDiffProteome)

                                        pred.append(result)

                for f in concurrent.futures.as_completed(pred):

                        epitome.append(f.result())


        if ic50cutoff:

                epitome = [element for subArray in epitome for element in subArray if float(element[3]) < ic50cutoff]  # discard epitopes with high predicted IC50

        else:

                epitome = [element for subArray in epitome for element in subArray]

        [element.append(predictImmu([str(element[2])], str(element[1]))[2]) for element in epitome]

        [element.insert(0, transcriptsToGeneMap[element[0]]) for element in epitome]


        return epitome



def calcGPIE(data):

        data[4] = pd.to_numeric(data[4], downcast = "float")
        data[5] = pd.to_numeric(data[5], downcast = "float")
        data['f4'] = pd.to_numeric(data['f4'], downcast = "float")

        f1 = (data[4] - data[4].max()) / (data[4].min() - data[4].max())
        f2 = (data[5] - data[5].min()) / (data[5].max() - data[5].min())
        f3 = data['TumorMedian'] / 100
        f3[f3 > 1] = 1

        data['gPie'] = 100 * f1 * f2* f3 * data['f4']

        return data

def calcImmunochem(input_pep):
	
	analysed_seq = ProteinAnalysis(input_pep)
	molW = analysed_seq.molecular_weight()
	stability = analysed_seq.instability_index()
	iso = analysed_seq.isoelectric_point()
	hydrophobicity = analysed_seq.gravy()

	aa_polarity = {'A':0.00, 'R':52.00, 'D':49.70, 'N':3.38, 'C': 1.48,'E':49.90, 'Q':3.53, 'G':0.00, 'H':51.60, 'L':0.13, 'I' : 0.13, 'K':49.50, 'M':1.43, 'F':0.35, 'P':1.58, 'S':1.67, 'T':1.66, 'W':2.10,'Y':1.61, 'V':0.13}

	polarity = 0

	for amino_acid in input_pep:
		if amino_acid in aa_polarity:

			polarity += aa_polarity[amino_acid]

	polarity /= len(input_pep)
	
	physChem = [input_pep, iso, molW, stability, hydrophobicity, polarity]

	return physChem
