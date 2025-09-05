import shutil
from UTILS import *

def convertVCFToBfile(VCF, name, folder, plink1, logFile, debug):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --vcf {VCF} --make-bed --out {outputPrefix} --keep-allele-order --double-id"

    execute(commandLine, logFile, debug)
    saveInfoLog(f"{outputPrefix}", "plink1", "ConvertVCFToPlink", logFile)

    return outputPrefix

def extractSamples(bfile, extractFile, name, folder, plink1, logFile, debug):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {bfile} --make-bed --out {outputPrefix} --keep {extractFile}"

    execute(commandLine, logFile, debug)
    saveInfoLog(f"{outputPrefix}", "plink1", "ExtractSamples", logFile)

    return outputPrefix

def extractVariants(bfile, extractFile, name, folder, plink1, logFile, debug):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {bfile} --make-bed --out {outputPrefix} --extract {extractFile}"
    execute(commandLine,  logFile, debug)
    saveInfoLog(f"{outputPrefix}", "plink1", "ExtractVariants", logFile)

    return outputPrefix

def mergeCommonFiles(file1, file2, name, folder, plink1, logFile, debug):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {file1} --bmerge {file2} --make-bed --out {outputPrefix}"
    execute(commandLine, logFile, debug)
    saveInfoLog(f"{outputPrefix}", "plink1", "MergeRefAndAlt", logFile)

    return outputPrefix
def mergeRefAndTarget(bfileTarget, bfileRef, folder, name, plink1, logFile, debug):
    print(f"Merging reference and target")

    createFolder(folder, False)

    command = f"{plink1} --bfile {bfileTarget} --make-bed --out {folder}/TargetMAF --maf 0.01"
    execute(command, logFile, debug)
    saveInfoLog(f"{folder}/TargetMAF", "plink1", "MAFTarget", logFile)


    shutil.copy2(f"{folder}/TargetMAF.bed",f"{folder}/TargetMAFID.bed")
    shutil.copy2(f"{folder}/TargetMAF.fam", f"{folder}/TargetMAFID.fam")

    command = f"{plink1} --bfile {bfileRef} --make-bed --out {folder}/RefMAF --maf 0.01"
    execute(command, logFile, debug)
    saveInfoLog(f"{folder}/TargetMAF", "plink1", "MAFRef", logFile)

    shutil.copy2(f"{folder}/RefMAF.bed", f"{folder}/RefMAFID.bed")
    shutil.copy2(f"{folder}/RefMAF.fam", f"{folder}/RefMAFID.fam")


    print(f"Changing the bim ID for target file")
    #Open bim file to change varID to chrom:pos:A1:A2
    dictVar = {}
    bimFileIn = open(f"{folder}/TargetMAF.bim")
    bimFileOut = open(f"{folder}/TargetMAFID.bim", "w")
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
    bimFileIn = open(f"{folder}/RefMAF.bim")

    bimFileOut = open(f"{folder}/RefMAFID.bim", "w")
    inCommon = open(f"{folder}/inCommonRefAlt.txt", "w")
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

    refFile = f"{folder}/RefMAFID"
    targetFile = f"{folder}/TargetMAFID"
    commonList = f"{folder}/inCommonRefAlt.txt"
    targetCommon = extractVariants(targetFile, commonList, f"{name}_Target_common", f"{folder}/", plink1, logFile, debug)
    refCommon = extractVariants(refFile, commonList, f"{name}_Ref_common", f"{folder}/", plink1, logFile, debug)

    return mergeCommonFiles(targetCommon, refCommon, f"{name}_Merged", f"{folder}/", plink1, logFile, debug), targetCommon, refCommon


def separateGenotypedDataBySex(genotyped, covarDict, folder, remove, name, plink1, logFile, debug):
    print("Split the genotyping data based on the declared sex on the covar file\n")

    bfile = genotyped
    if genotyped.endswith("vcf") or genotyped.endswith("vcf.gz"):
        bfile = convertVCFToBfile(genotyped, f'toExtractSex', folder, plink1, logFile, debug)

    extractFileFemale = f"{folder}/{name}_Female_toExtract.txt"
    extractFileMale = f"{folder}/{name}_Male_toExtract.txt"

    print(f"\tCreating the files to split sexes ({extractFileMale} and {extractFileFemale})")
    fileMale = open(extractFileMale, "w")
    fileFemale = open(extractFileFemale, "w")

    for sample in covarDict:
        if covarDict[sample]["SEX"] == "1":
            fileMale.write(f"{sample} {sample}\n")
        if covarDict[sample]["SEX"] == "2":
            fileFemale.write(f"{sample} {sample}\n")


    fileMale.close()
    fileFemale.close()

    if remove:
        fileRemove = open(f"{folder}/outliersListMaleFemale.txt", "w")
        for ID in remove:
            fileRemove.write(f"{ID}\t{ID}\n")
        fileRemove.close()

        bfile = removeSamples(bfile, f"{folder}/outliersListMaleFemale.txt", f"{name}_Both", folder, plink1, logFile, debug)


    bfileFemale = extractSamples(bfile, extractFileFemale, f"{name}_Female", folder, plink1, logFile, debug)
    bfileMale = extractSamples(bfile, extractFileMale, f"{name}_Male", folder, plink1, logFile, debug)

    return bfileFemale, bfileMale, bfile

def removeSamples(bfile, removeFile, name, folder, plink1, logFile, debug, run= True):
    outputPrefix = f"{folder}/{name}"
    commandLine = f"{plink1} --bfile {bfile} --make-bed --out {outputPrefix} --remove {removeFile}"
    execute(commandLine,  logFile, debug, run)
    saveInfoLog(f"{outputPrefix}", "plink1", "RemoveVariants", logFile)

    return outputPrefix
def removeLDAndMAFToPCARefTarget(merged, target, ref, folder, name, plink1, logFile, debug):
    command = f"{plink1} --bfile {merged} --indep-pairwise 200 50 0.2 --out {folder}/{name} --maf 0.01"
    execute(command, logFile, debug)

    targetLD = extractVariants(target, f"{folder}/{name}.prune.in", f"{name}_Target", folder, plink1, logFile, debug)
    refLD = extractVariants(ref, f"{folder}/{name}.prune.in", f"{name}_Ref", folder, plink1, logFile, debug)

    return targetLD, refLD

def removeLDAndMAFToPCARef(target, folder, name, plink1, logFile, debug):
    command = f"{plink1} --bfile {target} --indep-pairwise 200 50 0.2 --out {folder}/{name} --maf 0.01"
    execute(command, logFile, debug)

    return extractVariants(target, f"{folder}/{name}.prune.in", f"{name}_Target", folder, plink1, logFile, debug)


def imputeSex(bfile, name, folder, plink1, logFile, debug):
    outputPrefix = f"{folder}/{name}"

    commandLine = (f"{plink1} --bfile {bfile} --make-bed --out {outputPrefix}_Split --split-x hg38 no-fail")
    execute(commandLine, logFile, debug)
    saveInfoLog(f"{outputPrefix}_Split", "plink1", "Split X (ImputeSexRef)", logFile)

    commandLine = (f"{plink1} --bfile {outputPrefix}_Split --make-bed --out {outputPrefix} --impute-sex 0.4 0.8")
    execute(commandLine, logFile, debug)
    saveInfoLog(f"{outputPrefix}_Split", "plink1", "Impute-Sex (ImputeSexRef)", logFile)

    return outputPrefix

def filterAutosomal(X, autosomal, name, folder, plink, logFile, debug):
    listKeep = open(f"{folder}/toExtractFromAutosomal.txt", "w")
    fam = open(f"{X}.fam")

    for line in fam:
        FID, IID, M, F, Se, St = line.strip().split()
        listKeep.write(f"{IID}\t{IID}\n")
    listKeep.close()

    return extractSamples(autosomal, f"{folder}/toExtractFromAutosomal.txt", name, folder, plink, logFile, debug)


def convertVCFToPFile(covarFile, config, folder, name, update, logFile, debug):
    imputedVCF = config["input"]["x"]["imputed"]
    plink2 = config["programs"]["plink2"]

    fileIn = open(covarFile)
    fileExtract = open(f"{folder}/{name}_toExtract.txt", "w")
    header = True
    for line in fileIn:
        if header:
            header = False
        else:
            ID = line.split()[0]
            fileExtract.write(f"{ID}\n")
    fileExtract.close()

    line = f"{plink2} --vcf {imputedVCF} --update-sex {update} --keep {folder}/{name}_toExtract.txt --make-pgen --out {folder}/{name}"
    execute(line, logFile, debug)
    return f"{folder}/{name}"
