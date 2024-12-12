import sys
import numpy as np
sys.path.append('../')
from Handlers.PLINK1 import *


def createOutlierList(covarDict, bfile, PCASource, PCList, folder, name):
    print(f"{bfile}.fam")
    print(f"{covarDict}")
    famFile = open(f"{bfile}.fam")
    sampleList = []

    for line in famFile:
        FID, IID, mother, father, sex, pheno = line.strip().split()
        sampleList.append(IID)

    removalDict = {}

    for dataSource in PCASource:
        for PC in PCList:
            listPC = []
            for ind in sampleList:
                if ind in covarDict:
                    PCInfo = float(covarDict[ind]["PCA"][dataSource][PC])
                    listPC.append(PCInfo)

            mean = np.mean(listPC)
            sd = np.std(listPC)

            lowerBound = mean - 3 * sd
            upperBound = mean + 3 * sd

            for ind in sampleList:
                if ind in covarDict:
                    PCInfo = float(covarDict[ind]["PCA"][dataSource][PC])
                    if PCInfo < lowerBound or PCInfo > upperBound:
                        if ind not in removalDict:
                            removalDict[ind] = []
                        removalDict[ind].append(f"{PC}({dataSource})")

    fileToRemove = open(f"{folder}/{name}_outlierList.txt", "w")
    fileToRemoveExplanation = open(f"{folder}/{name}_outlierExplanation.txt", "w")

    for ind in removalDict:
        fileToRemove.write(f"{ind}\t{ind}\n")
        fileToRemoveExplanation.write(f"{ind}:")
        for info in removalDict[ind]:
            fileToRemoveExplanation.write(f" {info}")
        fileToRemoveExplanation.write(f"\n")

    fileToRemove.close()
    fileToRemoveExplanation.close()

    return f"{folder}/{name}_outlierList.txt"

def getProjectedPCA(target, reference, gcta, X, threads, name, folder, plink1, logFile):
    command = f"{gcta} --bfile {reference} --maf 0.01 --make-grm --out {folder}/{name} --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile)

    command = f"{gcta} --grm {folder}/{name} --maf 0.01 --pca 50 --out {folder}/{name}_PCs --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile)

    command = f"{gcta} --bfile {reference} --maf 0.01 --pc-loading {folder}/{name}_PCs --out {folder}/{name}_VariantLoading --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile)

    command = f"{gcta} --bfile {target} --maf 0.01 --project-loading {folder}/{name}_VariantLoading 50 --out {folder}/{name}_TargetPC --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile)

    return f"{folder}/{name}_TargetPC.proj.eigenvec"

def getPCA(bfile, PCA, name, folder, plink1, logFile):
    fileNonLD = removeLD(bfile, name, folder, plink1, logFile)

    command = f"Rscript {PCA} {fileNonLD} {folder}/{name}.gds {folder}/{name}_PCA.tsv"
    execute(command, logFile)

    return f"{folder}/{name}_PCA.tsv"

def execute(commandLine, logFile):
    os.system(commandLine)

    writeInfoFile = True
    logFile.write("Command line:\n")
    logFile.write(f"\t{commandLine}:\n")
    logFile.write(f"\n")
    logFile.write("=================================================================================================\n")
