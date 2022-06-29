#!/usr/bin/env Rscript

# args[1] -> Path to epitope predictions
# args[2] -> Path to file with network interaction scores
# args[3] -> Path to out path

library(randomForest)
library(scales)  # for color manipulation with alpha()


args = commandArgs(trailingOnly = T)

# Check if arguments supplied
if (length(args) < 3) {
  stop("Not enough arguments found!", call.=FALSE)
} else {
  epi.file = args[1]
  nis.file = args[2]
  out.path = args[3]
  if (file.access(epi.file, 4) != 0)
    stop(paste("Error: Cannot read file", epi.file))
  else
    cat("Found epitope predictions: ", epi.file, "\n")
  if (file.access(nis.file, 4) != 0)
    stop(paste("Error: Cannot read file", nis.file))
  else
    cat("Found network interaction scores: ", nis.file, "\n")
  if (file.access(out.path, 2) != 0)
    stop(paste("Error: Cannot write into", out.path))
  else
    cat("Preparing outputs in: ", out.path, "\n")
}


# New data, output from the curatopes pipeline
# Replace the path to input from *ARGS*
Results = read.table(gzfile(epi.file))

# add interaction score from network analysis
nis = read.table(nis.file, sep = "\t", header = T)  # network interaction score
Results$NetworkInteraction = nis[match(Results$GeneID, nis$ENSEMBL), "normalized"]

# The models were trained on features inherent to the peptide (no HLA etc),
# so we predict once per unique peptide.
Peptides = Results[!duplicated(Results[, "Peptide"]), ]
Peptides$Length = nchar(Peptides$Peptide)

# Activity
ActivityModels = readRDS(file="../Data/Enviroments/Balanced-Activity-Models.rds")


for(i in seq_along(BindingModels)){
  
  cat(paste("Processing Binding Model", i, "\n"))
  
  name = paste0("Model", i)
  Peptides[name] = as.numeric(predict(BindingModels[[i]], Peptides)) - 1

}

Peptides$BindingProbability =  rowMeans(Peptides[, grep("^Model[0-9]+$", colnames(Peptides))])

Peptides = Peptides[, -grep("^Model[0-9]+$", colnames(Peptides))]


# Binding
BindingModels = readRDS(file="../Data/Enviroments/Weighted-Binding-Models.rds")


for(i in seq_along(ActivityModels)){
  
  cat(paste("Processing Activity Model", i, "\n"))
  
  name = paste0("Model", i)
  Peptides[name] = as.numeric(predict(ActivityModels[[i]], Peptides)) - 1
}

Peptides$ActivityProbability =  rowMeans(Peptides[, grep("^Model[0-9]+$", colnames(Peptides))])

Peptides = Peptides[, -grep("^Model[0-9]+$", colnames(Peptides))]

# plot model output against each other
png(file.path(out.path, "BindingVsActivity_unique-peptides.png"), width = 960, height = 960)
plot(ActivityProbability ~ BindingProbability, Peptides, col = alpha("red", .05), pch = 16, cex = 2, xlim = c(0, 1) , ylim = c(0, 1))
dev.off()

Results = merge(Results, Peptides[, c("Peptide", "Length", "BindingProbability", "ActivityProbability")], by = "Peptide")


# helper function for score calculations
linear.norm = function(x, mi = min(x), ma = max(x), lower.cap = NA, upper.cap = NA)  {
  x = (x - mi)/(ma - mi)
  x[x < lower.cap] = lower.cap
  x[x > upper.cap] = upper.cap
  x
}

# calculate gPIE
Results$Scorev1 = with(Results, {(
  100
  * linear.norm(-IC50)
  * linear.norm(Immunogenicity)
  * linear.norm(TumorMedian, mi = 0, ma = 100, upper.cap = 1)
  * ExpressionIndex
)})

# calculate new score
Results$Scorev2 = with(Results, {(
  100
  * linear.norm(-log10(IC50), mi = -log10(2000), ma = -log10(30), lower.cap = 0, upper.cap = 1)
  * BindingProbability
  * ActivityProbability
  * NetworkInteraction
  * linear.norm(TumorMedian, mi = 0, ma = 100, upper.cap = 1)
  * ExpressionIndex
)})

Results$IC50.log10 = log10(Results$IC50)
Results$TumorMedian.log10 = log10(Results$TumorMedian)
Trifecta = Results[!duplicated(Results[, c("GeneID", "Peptide", "HLA.Allele")]), ]
Epitopes = Trifecta[!duplicated(Trifecta[, c("Peptide", "HLA.Allele")]), ]

plot.parameters.peptide = c("Length", "MolW", "Iso", "Stability", "Hydrophobicity", "Polarity", "BindingProbability", "ActivityProbability")
plot.parameters.epitope = c("IC50.log10")
plot.parameters.gene = c("NetworkInteraction", "TumorMedian.log10", "ExpressionIndex")
plot.parameters.scores = c("Scorev1", "Scorev2")
# densitiy plots on peptide level
for (pp in c(plot.parameters.peptide)) {
  png(file.path(out.path, paste0(pp, "-density_peptide-level.png")), width = 960, height = 960)
  plot(density(Peptides[[pp]]), main = pp)
  dev.off()
}
# density plots on epitope level
for (pp in c(plot.parameters.epitope)) {
  png(file.path(out.path, paste0(pp, "-density_epitope-level.png")), width = 960, height = 960)
  plot(density(Epitopes[[pp]]), main = pp)
  dev.off()
}
# density plots on gene-epitope level
for (pp in c(plot.parameters.gene, plot.parameters.scores)) {
  png(file.path(out.path, paste0(pp, "-density_gene-level.png")), width = 960, height = 960)
  plot(density(Trifecta[[pp]]), main = pp)
  dev.off()

  png(file.path(out.path, paste0(pp, "-density_gene-level_disregard-0.png")), width = 960, height = 960)
  plot(density(Trifecta[[pp]][Trifecta[[pp]] != 0]), main = pp)
  dev.off()
}
# plot against new score
plot.parameters = c(plot.parameters.peptide, plot.parameters.epitope, plot.parameters.gene)
for (pp in plot.parameters) {
  png(file.path(out.path, paste0("ScoreV2vs", pp, "_gene-level.png")), width = 960, height = 960)
  plot(as.formula(paste("Scorev2 ~", pp)), Trifecta, col = alpha("red", .05), pch = 16, cex = 2, ylim = c(0, 100), ylab = "ScoreV2")
  dev.off()
}
# plot new score against old score
png(file.path(out.path, "ScoreV2vsScoreV1.png"), width = 960, height = 960)
plot(Scorev2 ~ Scorev1, Trifecta, col = alpha("red", .05), pch = 16, cex = 2, xlim = c(0, 100), ylim = c(0, 100), xlab = "ScoreV1 (gPIE)", ylab = "ScoreV2")
dev.off()
# plot historgram of new score where it counts
png("ScoreV2-hist_greater10.png", width = 960, height = 960)
hist(subset(Trifecta, Scorev2 > 10)$Scorev2, breaks = 10:ceiling(max(Trifecta$Scorev2)), main = "Candidates with score > 10", cex.main = 2, cex.axis = 2)
dev.off()

write.table(Results, file = gzfile(file.path(out.path, "Results-scored.csv.gz")), sep = "\t", quote = F, row.names = F)
write.table(
  Trifecta[, !colnames(Trifecta) %in% c("TranscriptID", "IC50.log10", "TumorMedian.log10")],
  file = gzfile(file.path(out.path, "Results-scored-unique.csv.gz")), sep = "\t", quote = F, row.names = F
)
cat("Scores computed.\n")

