import sys
sys.path.append('../')
from Utils.GENERAL import execute

def countFam(filePrefix):
    famFile = open(f"{filePrefix}.fam")
    count = 0

    for line in famFile:
        count = count+1
    famFile.close()
    return count

def countBim(filePrefix):
    bimFile = open(f"{filePrefix}.bim")
    count = 0

    for line in bimFile:
        count = count+1
    bimFile.close()
    return count

def convertVCFToBfile(VCF, name, folder, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --vcf {VCF} --make-bed --out {outputPrefix} --double-id --keep-allele-order"

    execute(commandLine, "bfile", outputPrefix, logFile)

    return outputPrefix

def extractSamples(bfile, extractFile, name, folder, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {bfile} --make-bed --out {outputPrefix} --extract {extractFile}"

    execute(commandLine, "bfile", outputPrefix, logFile)

    return outputPrefix

def imputeSex(bfile, name, folder, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {bfile} --make-bed --out {outputPrefix} --impute-sex 0.4 0.8"

    execute(commandLine, "bfile", outputPrefix, logFile)

    return outputPrefix

def separateGenotypedDataBySex(genotyped, covarDict, name, folder, plink1, logFile):
    print("Split the genotyping data based on the declared sex on the covar file\n")

    bfile = convertVCFToBfile(genotyped, f'{name}_toExtractSex', folder, plink1, logFile)

    extractFileFemale = f"{folder}/{name}_Female_toExtract.txt"
    extractFileMale = f"{folder}/{name}_Male_toExtract.txt"

    fileMale = open(extractFileMale, "w")
    fileFemale = open(extractFileFemale, "w")

    for sample in covarDict:
        if covarDict[sample]["SEX"] == "1":
            fileMale.write(f"{sample} {sample}")
        if covarDict[sample]["SEX"] == "2":
            fileFemale.write(f"{sample} {sample}")

    fileMale.close()
    fileFemale.close()

    bfileFemale = extractSamples(bfile, extractFileFemale, f"{name}_Female", folder, plink1, logFile)
    bfileMale = extractSamples(bfile, extractFileMale, f"{name}_Male", folder, plink1, logFile)

    return bfileFemale, bfileMale

def separateGenotypedDataBySexWithoutCovar(genotyped, name, folder, plink1, logFile):
    print("To reference file to PCA we will infer the sex using PLINK --impute-sex\n")

    bfile = convertVCFToBfile(genotyped, f'{name}_toExtractSex', folder, plink1, logFile)
    bfileSex = imputeSex(bfile, f"{name}_imputeSex", folder, plink1, logFile)

    extractFileFemale = f"{folder}/{name}_Female_toExtract.txt"
    extractFileMale = f"{folder}/{name}_Male_toExtract.txt"

    fileMale = open(extractFileMale, "w")
    fileFemale = open(extractFileFemale, "w")

    famFile = open(f"{bfileSex}.fam")

    for line in famFile:
        FID, IID, mother, father, sex, pheno  = line.strip().split()

        if sex == "1" or sex == 1:
            fileMale.write(f"{IID}\t{FID}\n")
        elif sex == "2" or sex == 2:
            fileFemale.write(f"{IID}\t{FID}\n")
        else:
            print(f"PLINK was not able to infer the sex of  {IID} (cutoff: female (max) 0.6, male (min) 0.8)\n")

    fileMale.close()
    fileFemale.close()

    bfileFemale = extractSamples(bfile, extractFileFemale, f"{name}_Female", folder, plink1, logFile)
    bfileMale = extractSamples(bfile, extractFileMale, f"{name}_Male", folder, plink1, logFile)

    return bfileFemale, bfileMale

def mergeRefAndTarget(bfileTarget, bfileRef, name, plink1, logFile):
    pass