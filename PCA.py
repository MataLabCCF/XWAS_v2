import numpy as np
from UTILS import *
from PLINK import *

def getNonProjectedPCA(target, gcta, X, threads, name, folder, logFile, debug):
    # Projected
    command = f"{gcta} --bfile {target} --maf 0.01 --make-grm --out {folder}/{name} --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile, debug)

    command = f"{gcta} --grm {folder}/{name} --maf 0.01 --pca 50 --out {folder}/{name}_PCs --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile, debug)
    return f"{folder}/{name}_PCs.eigenvec"

def getProjectedPCA(target, reference, gcta, X, threads, name, folder, logFile, debug):
    print("Calculating the projected PCA")
    command = f"{gcta} --bfile {reference} --maf 0.01 --make-grm --out {folder}/{name} --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile, debug)

    command = f"{gcta} --grm {folder}/{name} --maf 0.01 --pca 50 --out {folder}/{name}_PCs --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile, debug)

    command = f"{gcta} --bfile {reference} --maf 0.01 --pc-loading {folder}/{name}_PCs --out {folder}/{name}_VariantLoading --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile, debug)

    command = f"{gcta} --bfile {target} --maf 0.01 --project-loading {folder}/{name}_VariantLoading 50 --out {folder}/{name}_TargetPC --thread-num {threads}"
    if X:
        command = f"{command} --chr 23 --autosome-num 25"
    else:
        command = f"{command} --autosome"
    execute(command, logFile, debug)

    print(f"Returning the file: {folder}/{name}_TargetPC.proj.eigenvec")
    return f"{folder}/{name}_TargetPC.proj.eigenvec"


def LDAndPCA(bfile, name, folder, plink, gcta, X, threads, logFile, debug):
    targetLD = removeLDAndMAFToPCARef(bfile, folder, name, plink, logFile, debug)

    return getNonProjectedPCA(targetLD, gcta, X, threads, name, folder, logFile, debug)

def mergeAndLDAndProjectPCA(target, refAll, extract, name, folder, plink, gcta, X, threads, logFile, debug):
    refSex = refAll
    if extract:
        refSex = extractSamples(refAll, extract, f"Ref_{name}", folder, plink, logFile, debug)
    mergedCommon, targetCommon, refCommon = mergeRefAndTarget(target, refSex,f"{folder}/Merge", f"{name}ToProject", plink, logFile, debug)
    targetLD, refLD = removeLDAndMAFToPCARefTarget(mergedCommon, targetCommon, refCommon, folder, name, plink, logFile, debug)

    return getProjectedPCA(targetLD, refLD, gcta, X, threads, name, folder, logFile, debug)

def getOutlier(fileName, outlierList, numPCs):
    print(f"\tOur outlier detection started with {len(outlierList)}. We will look for outliers on the first {numPCs} PCs on the file {fileName}")
    fileIn = open(fileName)
    dictPCA = {}
    sumStat = {}

    for line in fileIn:
        split = line.strip().split()
        for i in range(2, len(split)):
            PC = f"PC{i - 1}"
            if PC not in dictPCA:
                dictPCA[PC] = []
            dictPCA[PC].append(float(split[i]))
    fileIn.close()

    for PC in dictPCA:
        sumStat[PC] = {}
        sumStat[PC]["mean"] = np.mean(dictPCA[PC])
        sumStat[PC]["sd"] = np.std(dictPCA[PC])

    fileIn = open(fileName)
    for line in fileIn:
        split = line.strip().split()
        for i in range(2, len(split)):
            PCNum = i - 1
            if PCNum <= numPCs:
                PC = f"PC{PCNum}"
                PCInd = float(split[i])

                minRange = sumStat[PC]["mean"] - (3 * sumStat[PC]["sd"])
                maxRange = sumStat[PC]["mean"] + (3 * sumStat[PC]["sd"])

                if PCInd < minRange or PCInd > maxRange:
                    if split[1] not in outlierList:
                        outlierList.append(split[1])
    print(f"Leaving with {len(outlierList)} outliers on outlier list")
    return outlierList


def detectPCAOutliersMaleFemale(bfileMale, bfileFemale, config, folder, threads, logFile, debug):
    outlierList = []

    folderOutlierDetection = f"{folder}/OutlierDetection/"
    createFolder(folderOutlierDetection, False)

    plink = config["programs"]['plink1']
    gcta = config["programs"]['gcta']

    if "outlier" in config:
        print(f"Let's find the outliers")
        if "x" in config["outlier"]:
            if "Ref" in config["outlier"]["x"]:
                print(f"\t\tProjected X-PCA")
                if "reference" in config:
                    if "x" in config["reference"]:
                        bfileRef = config["reference"]['x']

                        #If VCF, convert to PLINK1
                        if bfileRef.endswith("vcf") or bfileRef.endswith("vcf.gz"):
                            bfileRef = convertVCFToBfile(bfileRef, f'referenceX', folderOutlierDetection, plink, logFile, debug)

                        # Impute sex to reference
                        refImputeSex = imputeSex(bfileRef, f"imputeSex", folderOutlierDetection, plink, logFile, debug)

                        # Make files to extract by sex
                        famFile = open(f"{refImputeSex}.fam")
                        extractMale = f"{folderOutlierDetection}/extractRefMale.txt"
                        extractFemale = f"{folderOutlierDetection}/extractRefFemale.txt"

                        fileMale = open(extractMale, "w")
                        fileFemale = open(extractFemale, "w")

                        for line in famFile:
                            FID, IID, Mo, Fa, sex, pheno = line.strip().split()
                            if sex == "1":
                               fileMale.write(f"{IID}\t{IID}\n")
                            elif sex == "2":
                                fileFemale.write(f"{IID}\t{IID}\n")
                            else:
                                print(f"{IID} had a F between 0.4 and 0.8, so not determined sex")

                        famFile.close()
                        fileMale.close()
                        fileFemale.close()

                        # Build in common files
                        PCAMale = mergeAndLDAndProjectPCA(bfileMale, bfileRef, extractMale, "ProjectedX-PCAMale", folderOutlierDetection, plink, gcta, True, threads, logFile, debug)
                        PCAFemale = mergeAndLDAndProjectPCA(bfileFemale, bfileRef, extractFemale, "ProjectedX-PCAFemale", folderOutlierDetection, plink, gcta, True, threads, logFile, debug)

                        print(f"Outlier Male {PCAMale} and outlier Female {PCAMale}")
                        outlierList = getOutlier(PCAMale, outlierList, config["outlier"]['pc']['x'])
                        outlierList = getOutlier(PCAFemale, outlierList, config["outlier"]['pc']['x'])
                    else:
                        sys.exit("Error: you requested Projected X-PCA to detect PCA outlier but you did not provide X reference")
                else:
                    sys.exit("Error: you requested Projected X-PCA to detect PCA outlier but you did not provide any reference (autosomal or X)")

            if "NonRef" in config["outlier"]["x"]:
                print(f"\t\tNon-projected X-PCA")

                PCAMale = LDAndPCA(bfileMale, "NonProjectedX-PCAMale", folderOutlierDetection, plink, gcta, True, threads, logFile, debug)
                PCAFemale = LDAndPCA(bfileFemale, "NonProjectedX-PCAFemale", folderOutlierDetection, plink, gcta, True, threads, logFile, debug)

                outlierList = getOutlier(PCAMale, outlierList, config["outlier"]['pc']['x'])
                outlierList = getOutlier(PCAFemale, outlierList, config["outlier"]['pc']['x'])
        if "autosomal" in config["outlier"]:
            if "Ref" in config["outlier"]["autosomal"]:
                print(f"\t\tProjected Autosomal PCA")
                if "reference" in config:
                    if "autosomal" in config["reference"]:
                        bfileRef = config["reference"]['autosomal']
                        bfileAuto = config["input"]["autosomal"]['typed']

                        #If VCF, convert to PLINK1
                        if bfileRef.endswith("vcf") or bfileRef.endswith("vcf.gz"):
                            bfileRef = convertVCFToBfile(bfileRef, f'referenceX', folderOutlierDetection, plink, logFile, debug)

                        extract = ""
                        PCAAutosomal = mergeAndLDAndProjectPCA(bfileAuto, bfileRef, extract, "ProjectedAutosomal", folderOutlierDetection, plink, gcta, False, threads, logFile, debug)
                        outlierList = getOutlier(PCAAutosomal, outlierList, config["outlier"]['pc']['autosomal'])
                    else:
                        sys.exit("Error: you requested Projected X-PCA to detect PCA outlier but you did not provide X reference")
                else:
                    sys.exit("Error: you requested Projected X-PCA to detect PCA outlier but you did not provide any reference (autosomal or X)")
            if "NonRef" in config["outlier"]["autosomal"]:
                print(f"\t\tNon-projected Autosomal")
                bfileAuto = config["input"]["autosomal"]['typed']

                PCAAutosomal = LDAndPCA(bfileAuto, "NonProjectedAutosomalPCA", folderOutlierDetection, plink, gcta, False, threads, logFile, debug)
                outlierList = getOutlier(PCAAutosomal, outlierList, config["outlier"]['pc']['x'])
    else:
        print(f"No outlier detection... Let's move foward")

    return outlierList

def getRegressionPCA(female, male, both, config, name, folder, threads, logFile, debug):
    plink = config["programs"]['plink1']
    gcta = config["programs"]['gcta']

    if "pca" in config:
        if "x" in config["pca"]:
            PCAFemale = LDAndPCA(female, "ToRegressionFemale", folder, plink, gcta, True, threads, logFile, debug)
            PCAMale = LDAndPCA(male, "ToRegressionMale", folder, plink, gcta, True, threads, logFile, debug)
            PCABoth = LDAndPCA(both, "ToRegressionBoth", folder, plink, gcta, True, threads, logFile, debug)

        else:
            if "autosomal" in config["pca"]:
                femaleAuto = filterAutosomal(female, config["input"]["autosomal"], f"{name}_female", f"{folder}/OutlierDetection/", plink, logFile, debug)
                maleAuto = filterAutosomal(male, config["input"]["autosomal"], f"{name}_male", f"{folder}/OutlierDetection/", plink, logFile, debug)
                bothAuto = filterAutosomal(both, config["input"]["autosomal"], f"{name}_both", f"{folder}/OutlierDetection/", plink, logFile, debug)

                PCAFemale = LDAndPCA(femaleAuto, "ToRegressionFemale", folder, plink, gcta, False, threads, logFile, debug)
                PCAMale = LDAndPCA(maleAuto, "ToRegressionMale", folder, plink, gcta, False, threads, logFile, debug)
                PCABoth = LDAndPCA(bothAuto, "ToRegressionBoth", folder, plink, gcta, False, threads, logFile, debug)


    else:
        sys.exit(f"We did not found the PCA flag in the config file")
    return PCAFemale, PCAMale, PCABoth

def detectPCAOutliersBoth(bfileBoth, config, folder, threads, logFile, debug):
    outlierList = []

    folderOutlierDetection = f"{folder}/OutlierDetection/"
    createFolder(folderOutlierDetection, False)

    plink = config["programs"]['plink1']
    gcta = config["programs"]['gcta']


    if "outlier" in config:
        print(f"Let's find the outliers")
        if "x" in config["outlier"]:
            if "Ref" in config["outlier"]["x"]:
                print(f"\t\tProjected X-PCA")
                if "reference" in config:
                    if "x" in config["reference"]:
                        bfileRef = config["reference"]['x']

                        #If VCF, convert to PLINK1
                        if bfileRef.endswith("vcf") or bfileRef.endswith("vcf.gz"):
                            bfileRef = convertVCFToBfile(bfileRef, f'referenceX', folderOutlierDetection, plink, logFile, debug)

                        PCABoth = mergeAndLDAndProjectPCA(bfileBoth, bfileRef, "","ProjectedX-PCABoth", folderOutlierDetection, plink, gcta, True, threads, logFile, debug)

                        print(f"Outlier Both {PCABoth}")
                        outlierList = getOutlier(PCABoth, outlierList, config["outlier"]['pc']['x'])

                    else:
                        sys.exit("Error: you requested Projected X-PCA to detect PCA outlier but you did not provide X reference")
                else:
                    sys.exit("Error: you requested Projected X-PCA to detect PCA outlier but you did not provide any reference (autosomal or X)")

            if "NonRef" in config["outlier"]["x"]:
                print(f"\t\tNon-projected X-PCA")

                PCABoth = LDAndPCA(bfileBoth, "NonProjectedX-PCABoth", folderOutlierDetection, plink, gcta, True, threads, logFile, debug)

                outlierList = getOutlier(PCABoth, outlierList, config["outlier"]['pc']['x'])
        if "autosomal" in config["outlier"]:
            # Building the list of samples from both for autosomal
            bfileAuto = config["input"]["autosomal"]['typed']

            fileFam = open(f"{bfileBoth}.fam")
            extract = open(f"{folder}/AllSamplesFromBoth.txt", "w")
            for line in fileFam:
                FID, IID, M, F, Sex, Pheno = line.strip().split()
                extract.write(f"{IID}\t{IID}\n")
            extract.close()
            bothAuto = extractSamples(bfileAuto, f"{folder}/AllSamplesFromBoth.txt", "BothSamples", f"{folder}/", plink, logFile, debug)

            if "Ref" in config["outlier"]["autosomal"]:
                print(f"\t\tProjected Autosomal PCA")
                if "reference" in config:
                    if "autosomal" in config["reference"]:
                        bfileRef = config["reference"]['autosomal']

                        #If VCF, convert to PLINK1
                        if bfileRef.endswith("vcf") or bfileRef.endswith("vcf.gz"):
                            bfileRef = convertVCFToBfile(bfileRef, f'referenceX', folderOutlierDetection, plink, logFile, debug)

                        extract = ""
                        PCAAutosomal = mergeAndLDAndProjectPCA(bothAuto, bfileRef, extract, "ProjectedAutosomalBoth", folderOutlierDetection, plink, gcta, False, threads, logFile, debug)
                        outlierList = getOutlier(PCAAutosomal, outlierList, config["outlier"]['pc']['autosomal'])
                    else:
                        sys.exit("Error: you requested Projected X-PCA to detect PCA outlier but you did not provide X reference")
                else:
                    sys.exit("Error: you requested Projected X-PCA to detect PCA outlier but you did not provide any reference (autosomal or X)")
            if "NonRef" in config["outlier"]["autosomal"]:
                print(f"\t\tNon-projected Autosomal")
                bfileAuto = config["input"]["autosomal"]['typed']

                PCAAutosomal = LDAndPCA(bothAuto, "NonProjectedAutosomalPCABoth", folderOutlierDetection, plink, gcta, False, threads, logFile, debug)
                outlierList = getOutlier(PCAAutosomal, outlierList, config["outlier"]['pc']['x'])
    else:
        print(f"No outlier detection... Let's move foward")

    return outlierList


