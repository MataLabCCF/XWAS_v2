import os
import gzip

def getUpdateSexFile(fileName, folder, name,covarDict):
    fileSex = open(f"{folder}/{name}_UpdateSex.txt", "w")
    fileSex.write(f"#IID\tSEX\n")

    vcfFile = gzip.open(fileName)

    header = True
    for line in vcfFile:
        line = line.decode("utf-8")
        if header:
            if "#CHROM" in line:
                split = line.strip().split()
                for i in range(9, len(split)):
                    if split[i] in covarDict:
                        if "SEX" in covarDict[split[i]]:
                            fileSex.write(f"{split[i]}\t{covarDict[split[i]]['SEX']}\n")
                        else:
                            fileSex.write(f"{split[i]}\t0\n")
                    else:
                        fileSex.write(f"{split[i]}\t0\n")

                break
    fileSex.close()
    return f"{folder}/{name}_UpdateSex.txt"



def convertVCFToPFile(VCF, bfile, outlier, covarDict, name, folder, plink2, logFile):
    #Geting the update sex file to PLINK2.
    #This update sex is the most annoying thing
    updateSexFile = getUpdateSexFile(VCF, folder, name, covarDict)

    outlierFile = open(outlier)
    outlierList = []

    #Get the outlier list to remove
    for line in outlierFile:
        FID, IID = line.strip().split()
        outlierList.append(IID)

    #Get the sample list without outliers to keep from imputed data
    famFile = open(f"{bfile}.fam")
    toKeepName = f"{folder}/{name}_toKeepToRegression"
    toKeepFile = open(toKeepName, "w")
    for line in famFile:
        FID, IID, mother, father, sex, pheno = line.strip().split()
        if IID not in outlierList:
            toKeepFile.write(f"{IID}\n")

    outlierFile.close()
    toKeepFile.close()

    #Conversion
    outputPrefix = f"{folder}/{name}"
    commandLine = (f"{plink2} --vcf {VCF} --make-pfile --out {outputPrefix} --keep {toKeepName} "
                   f"--update-sex {updateSexFile}")

    executePlink2(commandLine, "pfile", outputPrefix, logFile)

    os.system(f"mkdir {folder}/BackPFile")
    #os.system()


    return outputPrefix


def convertPfileToBfile(pfile, name, folder, plink2, logFile):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink2} --pfile {pfile} --make-bed --out {outputPrefix} --double-id --keep-allele-order"

    executePlink2(commandLine, "pfile", outputPrefix, logFile)

    return outputPrefix

def countPsam(filePrefix):
    famFile = open(f"{filePrefix}.psam")
    count = 0

    header = True
    for line in famFile:
        if header:
            header = False
        else:
            count = count + 1
    famFile.close()
    return count

def countPvar(filePrefix):
    bimFile = open(f"{filePrefix}.pvar")
    count = 0

    header = True
    for line in bimFile:
        if header:
            if "#CHROM" in line:
                header = False
        else:
            count = count + 1

    bimFile.close()
    return count


def executePlink2(commandLine, type, outputPrefix, logFile):
    os.system(commandLine)

    writeInfoFile = False

    if type == "pfile":
        writeInfoFile = True
        numSample = countPsam(outputPrefix)
        numVar = countPvar(outputPrefix)

    logFile.write("Command line:\n")
    logFile.write(f"\t{commandLine}:\n")
    logFile.write(f"\n")
    if writeInfoFile:
        logFile.write(f"Output information:\n")
        logFile.write(f"\tNumber of samples: {numSample}\n")
        logFile.write(f"\tNumber of varaintss: {numVar}\n")
        logFile.write(f"\n")
    logFile.write("=================================================================================================\n")