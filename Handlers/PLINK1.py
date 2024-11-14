import os
import sys

def removeLDAndMAFToPCA(merged, target, ref, folder, name, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    command = f"{plink1} --bfile {merged} --indep-pairwise 200 50 0.2 --out {folder}/{name} --maf 0.01"
    executePlink1(command, "in", outputPrefix, logFile)

    targetLD = extractVariants(target, f"{folder}/{name}.prune.in", f"{name}_Target", folder, plink1, logFile)
    refLD = extractVariants(ref, f"{folder}/{name}.prune.in", f"{name}_Ref", folder, plink1, logFile)

    return targetLD, refLD

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

    executePlink1(commandLine, "bfile", outputPrefix, logFile)

    return outputPrefix

def extractSamples(bfile, extractFile, name, folder, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {bfile} --make-bed --out {outputPrefix} --keep {extractFile}"

    executePlink1(commandLine, "bfile", outputPrefix, logFile)

    return outputPrefix

def removeLD(bfile, name, folder, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    command = f"{plink1} --bfile {bfile} --indep-pairwise 200 50 0.2 --out {folder}/{name}"
    executePlink1(command, "in", outputPrefix, logFile)
    return extractVariants(bfile, f"{folder}/{name}.prune.in", f"{name}_LD", folder, plink1, logFile)

def extractVariants(bfile, extractFile, name, folder, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {bfile} --make-bed --out {outputPrefix} --extract {extractFile}"

    executePlink1(commandLine, "bfile", outputPrefix, logFile)

    return outputPrefix

def mergeCommonFiles(file1, file2, name, folder, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {file1} --bmerge {file2} --make-bed --out {outputPrefix}"

    executePlink1(commandLine, "bfile", outputPrefix, logFile)
    return outputPrefix

def imputeSex(bfile, name, folder, plink1, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {bfile} --make-bed --out {outputPrefix} --impute-sex 0.4 0.8"

    executePlink1(commandLine, "bfile", outputPrefix, logFile)

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
            fileMale.write(f"{sample} {sample}\n")
        if covarDict[sample]["SEX"] == "2":
            fileFemale.write(f"{sample} {sample}\n")

    fileMale.close()
    fileFemale.close()

    bfileFemale = extractSamples(bfile, extractFileFemale, f"{name}_Female", folder, plink1, logFile)
    bfileMale = extractSamples(bfile, extractFileMale, f"{name}_Male", folder, plink1, logFile)

    return bfileFemale, bfileMale

def separateGenotypedDataBySexWithoutCovar(genotyped, name, folder, plink1, logFile):
    print("To reference file to PCA we will infer the sex using PLINK --impute-sex\n")

    if genotyped.endswith("vcf") or genotyped.endswith("vcf.gz"):
        bfile = convertVCFToBfile(genotyped, f'{name}_toExtractSex', folder, plink1, logFile)
    #else:
    #    bfile = convertPfileToBfile(genotyped, f'{name}_toExtractSex', folder, plink1, logFile)
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
            print(f"PLINK was not able to infer the sex of  {IID} (cutoff: female (max) 0.6, male (min) 0.8)")

    fileMale.close()
    fileFemale.close()

    bfileFemale = extractSamples(bfile, extractFileFemale, f"{name}_Female", folder, plink1, logFile)
    bfileMale = extractSamples(bfile, extractFileMale, f"{name}_Male", folder, plink1, logFile)

    return bfileFemale, bfileMale

def mergeRefAndTarget(bfileTarget, bfileRef, folder, name, plink1, logFile):
    print(f"Merging reference and target")

    os.system(f"mkdir {folder}/{name}_MergeRefAlt")

    command = f"{plink1} --bfile {bfileTarget} --make-bed --out {folder}/{name}_MergeRefAlt/TargetMAF --maf 0.01"
    executePlink1(command, "bfile", f"{folder}/{name}_MergeRefAlt/TargetMAF", logFile)
    os.system(f"cp {folder}/{name}_MergeRefAlt/TargetMAF.bed {folder}/{name}_MergeRefAlt/TargetMAFID.bed")
    os.system(f"cp {folder}/{name}_MergeRefAlt/TargetMAF.fam {folder}/{name}_MergeRefAlt/TargetMAFID.fam")


    os.system(f"cp {bfileRef}.bed {folder}/{name}_MergeRefAlt/Ref.bed")
    os.system(f"cp {bfileRef}.fam {folder}/{name}_MergeRefAlt/Ref.fam")


    print(f"Changing the bim ID for target file")
    #Open bim file to change varID to chrom:pos:A1:A2
    dictVar = {}
    bimFileIn = open(f"{folder}/{name}_MergeRefAlt/TargetMAF.bim")
    bimFileOut = open(f"{folder}/{name}_MergeRefAlt/TargetMAFID.bim", "w")
    for line in bimFileIn:
        chrom, ID, poscM, posBP, A1, A2 = line.strip().split()
        newID = f"{chrom}:{posBP}:{A1}:{A2}"
        bimFileOut.write(f"{chrom}\t{newID}\t{poscM}\t{posBP}\t{A1}\t{A2}\n")

        if newID not in dictVar:
            dictVar[newID] = True
        else:
            dictVar[newID] = False
    bimFileIn.close()
    bimFileOut.close()



    print(f"Changing the bim ID for ref file and generating the list of variants in common")
    # Open bim file to change varID to chrom:pos
    bimFileIn = open(f"{bfileRef}.bim")

    bimFileOut = open(f"{folder}/{name}_MergeRefAlt/Ref.bim", "w")
    inCommon = open(f"{folder}/{name}_MergeRefAlt/inCommon.txt", "w")
    for line in bimFileIn:
        chrom, ID, poscM, posBP, A1, A2 = line.strip().split()
        newID = f"{chrom}:{posBP}:{A1}:{A2}"
        bimFileOut.write(f"{chrom}\t{newID}\t{poscM}\t{posBP}\t{A1}\t{A2}\n")

        if newID in dictVar:
            #print(f"In common {newID}")
            if dictVar[newID]:
                #print(f"not added because it was false")
                inCommon.write(f"{newID}\n")

    inCommon.close()
    bimFileIn.close()
    bimFileOut.close()

    refFile = f"{folder}/{name}_MergeRefAlt/Ref"
    targetFile = f"{folder}/{name}_MergeRefAlt/TargetMAFID"
    commonList = f"{folder}/{name}_MergeRefAlt/inCommon.txt"
    targetCommon = extractVariants(targetFile, commonList, f"Target_common", f"{folder}/{name}_MergeRefAlt/", plink1, logFile)
    refCommon = extractVariants(refFile, commonList, f"Ref_common", f"{folder}/{name}_MergeRefAlt/", plink1, logFile)

    return mergeCommonFiles(targetCommon, refCommon, "Merged", f"{folder}/{name}_MergeRefAlt/", plink1, logFile), targetCommon, refCommon


def executePlink1(commandLine, type, outputPrefix, logFile):
    os.system(commandLine)

    writeInfoFile = False

    if type == "bfile":
        writeInfoFile = True
        numSample = countFam(outputPrefix)
        numVar = countBim(outputPrefix)

    logFile.write("Command line:\n")
    logFile.write(f"\t{commandLine}:\n")
    logFile.write(f"\n")
    if writeInfoFile:
        logFile.write(f"Output information:\n")
        logFile.write(f"\tNumber of samples: {numSample}\n")
        logFile.write(f"\tNumber of varaintss: {numVar}\n")
        logFile.write(f"\n")
    logFile.write("=================================================================================================\n")
