import os

def modelSelection(covar, modelSelection, folder, name, logFile):
    commandLine = f"Rscript {modelSelection} {covar} {folder}/{name}"

    covarFile = open(covar)
    for line in covarFile:
        allCovar = line.strip().split()
        break

    print("All covar in the covar file")
    print(allCovar)

    os.system(commandLine)
    logFile.write("Command line:\n")
    logFile.write(f"\t{commandLine}:\n")
    logFile.write(f"\n")
    logFile.write("=================================================================================================\n")

    print(f"Opening the file {folder}/{name}_variables.tsv")
    fileModel = open(f"{folder}/{name}_variables.tsv")

    model = []
    for line in fileModel:
        covarName = line.strip()
        if covarName in allCovar:
            print(f"Covar {covarName} being added")
            model.append(covarName)
        else:
            print(f"Covar {covarName} not present in the original covar list")
            for cov in allCovar:
                if cov in covarName:
                    print(f"\tCovar {covarName} is probably the changed version of {cov}")
                    if cov not in model:
                        print(f"\tCovar {cov} is added")
                        model.append(cov)

    print(f"Return the model: {model}")
    return model


def buildCovarFile(covarDict, outlier, PCList, PCASource, bfile, folder, name):
    outlierList = []

    outlierFile = open(outlier)
    for line in outlierFile:
        IID, FID = line.strip().split()
        outlierList.append(IID)
    outlierFile.close()

    covarList = []
    famFile = open(f"{bfile}.fam")
    covarFile = open(f"{folder}/{name}.tsv", "w")

    #Header -> IID <covar non PC> <Covar PC in Source_PCNumber format>
    covarFile.write("IID")
    for ind in covarDict:
        for covar in covarDict[ind]:
            if covar != "PCA":
                covarFile.write(f"\t{covar}")
                covarList.append(covar)
        break
    for source in PCASource:
        for PC in PCList:
            covarFile.write(f"\t{source}_{PC}")
    covarFile.write("\n")

    for line in famFile:
        FID, IID, mother, father, sex, pheno = line.strip().split()
        if IID not in outlierList and IID in covarDict:
            covarFile.write(f"{IID}")
            for covar in covarList:
                covarFile.write(f"\t{covarDict[IID][covar]}")
            for source in PCASource:
                for PC in PCList:
                    covarFile.write(f"\t{covarDict[IID]['PCA'][source][PC]}")
            covarFile.write("\n")
    covarFile.close()
    famFile.close()
    return f"{folder}/{name}.tsv"

def addPCAToCovarDict(filePCA, covarDict, dataSource):
    for ind in covarDict:
        #print(f"For error: {covarDict[ind]}")
        if "PCA" not in covarDict[ind]:
            covarDict[ind]["PCA"] = {}
        if dataSource not in covarDict[ind]["PCA"]:
            covarDict[ind]["PCA"][dataSource] = {}

    file = open(f"{filePCA}")

    #Project PCA has no header, while the previous PCA had
    PCList = []

    for line in file:
        data = line.strip().split("\t")
        ind = data[0]
        if ind in covarDict:
            for i in range(2, len(data)):
                PC = f"PC{i-1}"
                if PC not in PCList:
                    PCList.append(PC)
                covarDict[ind]["PCA"][dataSource][PC] = data[i]


    return covarDict, PCList



def readCovarFile(covarTable):
    print('We are reading the covar table. We are assuming that the ID is the first col')
    print('We are also assuming that there is the column SEX and the Phenotype column is named DISEASE')
    print('We are also checking if there is any covar field that is empty or NA')

    file = open(covarTable)

    covarDict = {}
    header = True

    # doNotFix is a flag to put the phenotype in PLINK2 format (1 control, 2 case)
    doNotFix = False
    for line in file:
        if header:
            header = False
            splitHeader = line.strip().split()
        else:
            split = line.strip().split()

            nonNA = True
            for i in range(len(split)):
                if split[i] == "NA" or split[i] == "" or split[i] == " " or split[i] == "nan":
                    nonNA = False
                    ind = split[0]
                    data = split[i]
                    headerName = splitHeader[i]

                    print(f'Removing the ind {ind} because there is missing data ({data}) on the field {headerName}')

            if nonNA:
                covarDict[split[0]] = {}
                for i in range(1, len(split)):
                    if splitHeader[i].upper() == "DISEASE":
                        if split[i] == "0" or split[i] == 0:
                            split[i] = "1"
                        elif split[i] == "1" or split[i] == 1:
                            split[i] = "2"
                        elif split[i] == "2" or split[i] == 2:
                            split[i] = "3"
                            doNotFix = True
                        else:
                            print(f"Unknown status value to sample {split[0]}: {split[i]} ")
                    covarDict[split[0]][splitHeader[i].upper()] = split[i]

    if doNotFix:
        for sample in covarDict:
            if covarDict[sample]["DISEASE"] == "2":
                covarDict[sample]["DISEASE"] = "1"
            elif covarDict[sample]["DISEASE"] == "3":
                covarDict[sample]["DISEASE"] = "2"
            else:
                print(f"Unknown sample DISEASE field:  {covarDict[sample]['DISEASE']}")

    return covarDict