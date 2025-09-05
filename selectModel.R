#Model selection
library(MASS)

options <- commandArgs(trailingOnly = TRUE)

csvName = options[1]
phenoName = options[2]
outputName = options[3]

#Reading other covar
info <- read.table(csvName, sep = "\t", header = T)

#Merging covar with PCA and removing the ID to variable selection
allCovarNonID = subset(info,select = -c(IID))

#Creating two models: Null and Full
formulaStart <- paste(phenoName, "~ .")
model <- glm(as.formula(formulaStart), data = allCovarNonID, family = "binomial")

#Selecting variables
X = stepAIC(model, direction = "both")
matrix = model.matrix(X)
names = colnames(matrix)[-1]

#Output
write.table(names, paste(outputName,"_variables.tsv", sep = ""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep = "\t")