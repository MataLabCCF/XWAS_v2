import sys
import numpy as np
sys.path.append('../')
from Handlers.PLINK1 import *

def createCovarFileToRegression(covarDict, filePCA, outFolder, fileName):
    print(f"Building covar file using {filePCA}")
    file = open(filePCA)
    fileOut = open(f"{outFolder}/{fileName}.tsv", "w")
    fileOut.write("IID")
    for ind in covarDict:
        for covar in covarDict[ind]:
            fileOut.write(f"\t{covar}")
        break
    for PC in range(1,51):
        fileOut.write(f"\tPC{PC}")
    fileOut.write(f"\n")


    for line in file:
        if not line.startswith("sample.id"):
            split = line.strip().split()
            if split[0] == split[1]:
                start = 2
            else:
                start = 1

            fileOut.write(f"{split[0]}")
            for field in covarDict[split[0]]:
                fileOut.write(f"\t{covarDict[split[0]][field]}")

            for i in range(start, len(split)):
                fileOut.write(f"\t{split[i]}")
            fileOut.write("\n")
    file.close()
    return f"{outFolder}/{fileName}.tsv"


def getOutlierFromPCA(filePCA, outlierList):
    print(f"Looking for outliers in {filePCA}")
    file = open(filePCA)

    dictPC = {}

    for line in file:
        if not line.startswith("sample.id"):
            split = line.strip().split()
            if split[0] == split[1]:
                start = 2
            else:
                start = 1

            j = 0
            for i in range(start, len(split)):
                j = j+1

                #Remove outlier in first 10 PCs
                if j <= 10:
                    if f"PC{j}" not in dictPC:
                        dictPC[f"PC{j}"] = {}
                        dictPC[f"PC{j}"]["dist"] = []
                    dictPC[f"PC{j}"]["dist"].append(float(split[i]))
    file.close()

    cutoff = {}
    for PC in dictPC:
        meanPC = np.mean(dictPC[PC]["dist"])
        sdPC = np.std(dictPC[PC]["dist"])

        cutoff[PC] = {}
        cutoff[PC]["mean"] = meanPC
        cutoff[PC]["min"] = meanPC-(3*sdPC)
        cutoff[PC]["max"] = meanPC+(3*sdPC)

    file = open(filePCA)

    for line in file:
        if not line.startswith("sample.id"):
            split = line.strip().split()
            if split[0] == split[1]:
                start = 2
            else:
                start = 1

            j = 0
            for i in range(start, len(split)):
                j = j + 1

                # Remove outlier in first 10 PCs
                if j <= 10:
                    minInterval = cutoff[f"PC{j}"]["min"]
                    maxInterval = cutoff[f"PC{j}"]["max"]
                    mean = cutoff[f"PC{j}"]["mean"]
                    if not (minInterval < float(split[i]) < maxInterval):
                        if split[0] not in outlierList:
                            print(f"{split[0]} removed PC{j} -> {minInterval} < {split[i]} < {maxInterval} MEAN: {mean} [{len(outlierList)}]")
                            outlierList.append(split[0])
    file.close()
    return outlierList


def getPCA(target, reference, gcta, X, threads, name, folder, plink1, logFile):
    print(f"\tCalculating PCA")
    #Non projectected
    if reference == "":
        # Projected
        command = f"{gcta} --bfile {target} --maf 0.01 --make-grm --out {folder}/{name} --thread-num {threads}"
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
        return f"{folder}/{name}_PCs.eigenvec"
    else:
        #Projected
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


def execute(commandLine, logFile):
    print(commandLine)
    os.system(commandLine)

    writeInfoFile = True
    logFile.write("Command line:\n")
    logFile.write(f"\t{commandLine}:\n")
    logFile.write(f"\n")
    logFile.write("=================================================================================================\n")
