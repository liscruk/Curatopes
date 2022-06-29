#!/usr/bin/env Rscript
args=commandArgs(TRUE)

### Logistic regression: Input file
TrainingSet = read.table("/home/lischecr/Projects/EpiPro/Scripts/TrainingSet.csv", sep = "\t",header=TRUE)
### Sampling
TrainingSet = TrainingSet[sample(1:nrow(TrainingSet)),]

### Sample Set aka Patient file
SampleSet = read.table(args[1],header = T, sep = "\t", na.strings = c(""))

#SampleSet=read.table("S437.ALL-ALLELES.Epitopes",header=T,sep = "\t")

Backup = SampleSet

### Negative scaling for features
TrainingSet$neg_IC50 = -(TrainingSet$IC50)
TrainingSet$neg_Polarity = -(TrainingSet$Polarity)
TrainingSet$neg_Mwt = -(TrainingSet$MolecularWeight)
TrainingSet$neg_Hydrophobicity = -(TrainingSet$Hydrophobicity)

### Extracting dataset for regression
data_set = subset(TrainingSet, select = c(2,4,6,7,12,13,14,15))

### Scaling and normalisation
norm_data_set = data_set
norm_data_set[,2:8] = scale(data_set[,2:8])

model = glm(Type ~., family=binomial(link = "logit"), data = norm_data_set)

### Fitting model on patient data
SampleSet = Backup

### Defining orders for features
SampleSet$Type = ("0")
SampleSet$neg_IC50 = -(SampleSet$IC50)
SampleSet$neg_Polarity = -(SampleSet$Polarity)
SampleSet$neg_Mwt = -(SampleSet$MolW)
SampleSet$neg_Hydrophobicity = -(SampleSet$Hydrophobicity)

### Extracting datset for regression
SampleSetModel = subset(SampleSet, select = c(19,13,11,18,20,21,22,23))
norm_SampleSet = SampleSetModel
norm_SampleSet$Immunogenicity=as.numeric(as.character(norm_SampleSet$Immunogenicity))
norm_SampleSet$Type = as.numeric(as.character(norm_SampleSet$Type))
norm_SampleSet[,2:8] = scale(norm_SampleSet[,2:8])

### Reloading the file for quick and dirty reset
SampleSet = Backup
SampleSet$glmPrediction = predict(model, newdata = subset(norm_SampleSet), type = 'response')

### Input variables
### X values
FPKM_1 = 1
FPKM_2 = 5
### Y values 
exp_1 = 0.1
exp_2 = 0.9
### Min/Max Normalisation vor VAF
SampleSet$VAF_norm = ((SampleSet$VAF - min(SampleSet$VAF))/(max(SampleSet$VAF)- min(SampleSet$VAF)))
### Equations
RatioFPKM = FPKM_1/FPKM_2
low_exp = (1-exp_1)/exp_1
high_exp = (1-exp_2)/exp_2
RatioExp = low_exp/high_exp
Inv_gm2_FPKM = -log10(RatioFPKM)/log10(RatioExp)
### Kinetic order
g_FPKM_VariantExp = 1/Inv_gm2_FPKM      
km_FPKM_VariantExp = FPKM_1*low_exp^Inv_gm2_FPKM
### Hill equation
SampleSet$score = SampleSet$VAF_norm * ((SampleSet$AVE^g_FPKM_VariantExp)/((km_FPKM_VariantExp^g_FPKM_VariantExp)+(SampleSet$AVE^g_FPKM_VariantExp)))
SampleSet$combinedScore = SampleSet$score*SampleSet$glmPrediction

Scored = SampleSet[order(-SampleSet$combinedScore),]
Scored = Scored[-c(13,14,15,18,19,20,21)]
colnames(Scored)[15] = c("Score")

write.table(Scored, "", sep = "\t", row.names = FALSE,quote = F)