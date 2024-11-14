library(GENESIS)
library(SeqArray)
library(SNPRelate)
library(GWASTools)


options <- commandArgs(trailingOnly = TRUE)

plinkName = options[1]
gdsName = options[2]
csvName = options[3]

bed = paste(plinkName,".bed", sep = "")
bim = paste(plinkName,".bim", sep = "")
fam = paste(plinkName,".fam", sep = "")

snpgdsBED2GDS(bed, fam, bim, gdsName)
genofile <- openfn.gds(gdsName)
pca<-snpgdsPCA(genofile,  autosome.only=FALSE)

tab <- data.frame(sample.id = pca$sample.id, PC1 = pca$eigenvect[,1], PC2 = pca$eigenvect[,2],
                  PC3 = pca$eigenvect[,3], PC4 = pca$eigenvect[,4], PC5 = pca$eigenvect[,5],
                  PC6 = pca$eigenvect[,6], PC7 = pca$eigenvect[,7], PC8 = pca$eigenvect[,8],
                  PC9 = pca$eigenvect[,9], PC10 = pca$eigenvect[,10], PC11 = pca$eigenvect[,11],
                  PC12 = pca$eigenvect[,12], PC13 = pca$eigenvect[,13], PC14 = pca$eigenvect[,14],
                  PC15 = pca$eigenvect[,15], PC16 = pca$eigenvect[,16], PC17 = pca$eigenvect[,17],
                  PC18 = pca$eigenvect[,18], PC19 = pca$eigenvect[,19], PC20 = pca$eigenvect[,20],
                  PC21 = pca$eigenvect[,21], PC22 = pca$eigenvect[,22], PC23 = pca$eigenvect[,23],
                  PC24 = pca$eigenvect[,24], PC25 = pca$eigenvect[,25], PC26 = pca$eigenvect[,26],
                  PC27 = pca$eigenvect[,27], PC28 = pca$eigenvect[,28], PC29 = pca$eigenvect[,29],
                  PC30 = pca$eigenvect[,30], PC31 = pca$eigenvect[,31], PC32 = pca$eigenvect[,32],
                  stringsAsFactors = FALSE)

write.table(tab, csvName, row.names=FALSE, quote=FALSE, sep = "\t")