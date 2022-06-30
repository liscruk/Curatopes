# <p align = center> Curatopes 1.5 </p>
## <p> A pipeline to generate self-tolerant, tumor restricted antigens from transcriptomics data </p>

Curatopes features a multi-level, multi-model prediction pipeline that allows, provided with tumor cohort transcriptomics data, to generate a list of possible tumor associciated antigens.

Currtenly, Curatopes is limited in its deployability since the pre-trained models are rather large.

Please contact us if you want to deploy the pipeline. We are working on making it more deployable.

General usage:
Add the `Curatopes.py` to your PATH to execute it from the command line.

Options:
```Python:
-e Provide a name for the tumor entity you are running (e.g Melanoma)
-i Provide an input file of tumor expression data in a gene by sample manner. 
   Expression should be in TPM or FPKM.   
   Format should be headerless with Ensemble IDs (Versionless) as Gene Identifiers.
-o Provide Path to outdir
-m Provide a list of MHC alleles of interest.
```

Dependencies can be supplied uppon request.

```Python:
#<<< <<< <<< Hardcoded files >>> >>> >>>#
# Use HPA low+medium+high
hpaList = "../Data/HPA/high+medium+low"
# GTEx tissue files need to be parsed separately
TissueStats = '../Data/GTEX/TissueStats'
# Proteome - ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
ProteomeAnnotation = "../Data/Proteome/Homo_sapiens.GRCh38.pep.all.protein_coding.txt"
# Critical Tissue List
CriticalTissue = "../Data/CriticalTissue/criticaltissues_binary.csv"
# network interaction scores
NIS = "../Data/Curatopes_UvealMelanoma/network-interaction-scores.csv"
```
