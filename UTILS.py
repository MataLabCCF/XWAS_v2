import os
import sys
import pandas as pd

def execute(line, logFile, debug, run = True):
    print("=======================================================================================================")
    print(line)
    if run:
        os.system(line)
    logFile.write(f"{line}\n")
    if debug:
        input(f"DEBUG MODE: PRESS ENTER")
    print("=======================================================================================================")

def saveInfoLog(inputName, type, analysis, logFile):
    if type == "plink1":
        inputFile = open(f"{inputName}.fam")

        count = 0
        for line in inputFile:
            count = count+1
        inputFile.close()
        logFile.write(f"Analysis {analysis}\n")
        logFile.write(f"\t#Samples {count}\n")

        inputFile = open(f"{inputName}.bim")
        count = 0
        for line in inputFile:
            count = count + 1
        inputFile.close()
        logFile.write(f"\t#Variants {count}\n")
        logFile.write(f"===============================================\n")



def readConfigFile(config):
    fileIn = open(config)

    dictConfig = {}
    covarDict = {}

    for line in fileIn:
        if not line.startswith("#"):
            fields = line.strip().split()

            if fields:

                if fields[0].lower() == "input":
                    if "input" not in dictConfig:
                        dictConfig["input"] = {}


                    if fields[1].lower() not in ["x", "autosomal", "covar"]:
                        sys.exit(f"Error: we do not recognize the {fields[1]} type on input field. We only recognize \'x\',"
                                 f" \'autosomal\' and \'covar\' (not case sensitive)")
                    else:
                        if fields[1].lower() not in dictConfig["input"]:
                            dictConfig["input"][fields[1].lower()] = {}
                        if fields[1].lower() == "x":
                            if fields[2].lower() == "typed":
                                dictConfig["input"][fields[1].lower()]["typed"] = fields[3]
                            elif fields[2].lower() == "imputed":
                                dictConfig["input"][fields[1].lower()]["imputed"] = fields[3]
                            else:
                                sys.exit(f"Error: we do not recognize the {fields[2]} type for X chromosome data. We just "
                                         f"accept \'typed\' (to PCA inferences) or \'imputed\' (to regression)")
                        elif fields[1].lower() == "autosomal":
                            if fields[2].lower() == "typed":
                                dictConfig["input"][fields[1].lower()]["typed"] = fields[3]
                            elif fields[2].lower() == "imputed":
                                sys.exit(f"Error: we do not recommend using imputed data to run PCA. If you want to use "
                                         f"imputed autosomal to PCA, please lie to us (replace the imputed to typed)")
                        elif fields[1].lower() == "covar":
                            print(f"Reading covar file ({fields[2]})")
                            covarDict = readCovarFile(fields[2], fields[3])
                            dictConfig["input"]["pheno"] = fields[3]

                if fields[0].lower() == "reference":
                    if "reference" not in dictConfig:
                        dictConfig["reference"] = {}

                    if fields[1].lower() == "x":
                        dictConfig["reference"][fields[1].lower()] = fields[2]
                    elif fields[1].lower() == "autosomal":
                        dictConfig["reference"][fields[1].lower()] = fields[2]
                    else:
                        sys.exit(f"Error: we do not recognize the {fields[2]} type for reference data. We just "
                                 f"accept \'X\' or \'autosomal\' (to regression)")

                if fields[0].lower() == "outlier":
                    if "outlier" not in dictConfig:
                        dictConfig["outlier"] = {}

                    if fields[1].lower() == "x":
                        if fields[1].lower() not in dictConfig["outlier"]:
                            dictConfig["outlier"][fields[1].lower()] = {}

                        if fields[2].lower() in ["y", "yes", "true"]:
                            dictConfig["outlier"][fields[1].lower()]["Ref"] = ""
                        elif fields[2].lower() in ["n", "no", "false"]:
                            dictConfig["outlier"][fields[1].lower()]["NonRef"] = ""
                        else:
                            sys.exit(f"Error: we do not recognize the {fields[2]} as an answer to use of projected X PCA. Please answer yes or no")

                    elif fields[1].lower() == "autosomal":
                        if fields[1].lower() not in dictConfig["outlier"]:
                            dictConfig["outlier"][fields[1].lower()] = {}

                        if fields[2].lower() in ["y", "yes", "true"]:
                            dictConfig["outlier"][fields[1].lower()]["Ref"] = ""
                        elif fields[2].lower() in ["n", "no", "false"]:
                            dictConfig["outlier"][fields[1].lower()]["NonRef"] = ""
                        else:
                            sys.exit(f"Error: we do not recognize the {fields[2]} as an answer to use of projected autosomal PCA. Please answer yes or no")
                    elif fields[1].lower() == "pc":
                        if fields[1].lower() not in dictConfig["outlier"]:
                            dictConfig["outlier"][fields[1].lower()] = {}

                        if fields[2].lower() == "x":
                            dictConfig["outlier"][fields[1].lower()]["x"] = int(fields[3])
                        elif fields[2].lower() == "autosomal":
                            dictConfig["outlier"][fields[1].lower()]["autosomal"] = int(fields[3])
                        else:
                            sys.exit(f"Error: we do not recognize {fields[2]} as an answer to define the number of PCs to be used on "
                                     f"outlier detection. Please answer with X or autosomal")
                    else:
                        sys.exit(f"Error: we do not recognize the {fields[1]} as an answer what PCA should be used to outlier detection. "
                                 f"We just accpet \'X\' or \'Autosomal\'")


                if fields[0].lower() == "pca":
                    if "pca" not in dictConfig:
                        dictConfig["pca"] = ""
                    else:
                        sys.exit(f"Error: we do not accept two PCA to regression. Please select X or autosomal")

                    if fields[1].lower() == "x":
                        dictConfig["pca"] = "x"
                    elif fields[1].lower() == "autosomal":
                        dictConfig["pca"] = "autosomal"
                    else:
                        sys.exit(
                            f"Error: we do not recognize the {fields[1]} as an answer what PCA should be used to regression. "
                            f"We just accpet \'X\' or \'Autosomal\'")

                if fields[0].lower() == "model":
                    if "model" not in dictConfig:
                        dictConfig["model"] = {}
                    else:
                        sys.exit(f"Error: we do not accept two models to regression. Please select between script ("
                                 f"and add the path to selectModel.R) or put your model (AGE+SEX+PC1)")

                    if fields[1].lower() == "script":
                        dictConfig["model"]["stepwise"] = fields[2]
                    elif fields[1].lower() == "regression":
                        dictConfig["model"]["finalmodel"] = fields[2]
                    else:
                        sys.exit("Error: Model has no valid option (script or regression)")

                if fields[0].lower() == "programs":
                    if "programs" not in dictConfig:
                        dictConfig["programs"] = {}

                    if fields[1].lower() == "plink1":
                        dictConfig["programs"]["plink1"] = fields[2]
                    elif fields[1].lower() == "plink2":
                        dictConfig["programs"]["plink2"] = fields[2]
                    elif fields[1].lower() == "gcta":
                        dictConfig["programs"]["gcta"] = fields[2]
                    elif fields[1].lower() == "gwama":
                        dictConfig["programs"]["gwama"] = fields[2]
                    elif fields[1].lower() == "r":
                        dictConfig["programs"]["r"] = fields[2]



    print(f"We finished to read the configuration file ")
    print(f"Your genetic inputs are: ")

    for data in dictConfig["input"]:
        if data != "covar" and data != "pheno":
            for type in dictConfig["input"][data]:
                print(f"\t{data}\t{dictConfig["input"][data][type]}({type})")
    if "reference" in dictConfig:
        print(f"We have references: ")
        for type in dictConfig["reference"]:
            print(f"\t{dictConfig["reference"][type]}({type})")
    if "outlier" in dictConfig:
        print(f"We will remove outliers: ")
        for type in dictConfig["outlier"]:
            if "Ref" in dictConfig["outlier"][type]:
                print(f"\tProjected {type} PCA... The reference is {dictConfig["reference"][type]}")
            if "NonRef" in dictConfig["outlier"][type]:
                print(f"\tRegular {type} PCA...")
    else:
        print(f"You are brave to run a XWAS without removing outliers. I liked that.")

    if "pca" in dictConfig:
        print(f"We will use {dictConfig['pca']} PCA on XWAS")
        if "stepwise" in dictConfig["model"]:
            print(f"\tThe PCs will be selected by the Stepwise regression implemented in the script: {dictConfig['model']}")
        else:
            print(f"\tThe model used will be {dictConfig['model']['finalmodel']}")

    else:
        print(f"Regression without PCA? Not on my house ")
        sys.exit("┏┓ \n"
                "┃┃╱╲ in\n"
                "┃╱╱╲╲ this\n"
                "╱╱╭╮╲╲house\n"
                "▔▏┗┛▕▔ we\n"
                "╱▔▔▔▔▔▔▔▔▔▔╲\n"
                "run XWAS with PCA\n"
                "╱╱┏┳┓╭╮┏┳┓ ╲╲\n"
                "▔▏┗┻┛┃┃┗┻┛▕▔\n")

    print(f"And all the tools are here:")
    for tool in ["plink1", "plink2", "gcta", "gwama", "r"]:
        print(f"\t{tool}: {dictConfig['programs'][tool]}")

    return dictConfig, covarDict
            

def createFolder(folder, log=True):
    if not os.path.exists(folder):
        os.mkdir(folder)
    if log:
        logFile = open(f"{folder}/XWAS.log", "w")
        return logFile



def readCovarFile(covarTable, phenoCol):
    print(f'We are reading the covar table. We are assuming that the ID is the first col and the phenotype is named {phenoCol}')
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
                    if splitHeader[i].upper() == phenoCol.upper():
                        if split[i] == "0" or split[i] == 0:
                            split[i] = "1"
                        elif split[i] == "1" or split[i] == 1:
                            split[i] = "2"
                        elif split[i] == "2" or split[i] == 2:
                            split[i] = "3"
                            doNotFix = True
                        else:
                            print(f"Unknown status value to sample {split[0]}: {split[i]} ")

                    if splitHeader[i].upper() == phenoCol.upper():
                        covarDict[split[0]]["DISEASE"] = split[i]
                    else:
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

def buildCovarFileAndModel(PCFile, covarDict, config, folder, name, logFile, debug):
    dictPCA = {}
    PCA = open(PCFile)
    for line in PCA:
        info = line.strip().split()
        ID = info[0]

        dictPCA[ID] = {}

        for i in range(2, len(info)):
            PC = f"PC{i - 1}"
            dictPCA[ID][PC] = info[i]
    PCA.close()

    fileOut = open(f"{folder}/{name}.tsv", "w")
    fileUpdate = open(f"{folder}/{name}_UpdateSex", "w")

    removalMissing = []
    outHeader = False

    for ID in dictPCA:
        if ID not in covarDict:
            print(f"Warning: {ID} not found on covar file")
        else:
            #Building header
            if not outHeader:
                fileOut.write(f"IID")
                for covar in covarDict[ID]:
                    fileOut.write(f"\t{covar.upper()}")
                for PC in range(1, i):
                    fileOut.write(f"\tPC{PC}")
                fileOut.write("\n")
                outHeader = True

            #Building file
            fileOut.write(f"{ID}")
            for covar in covarDict[ID]:
                if covar == "SEX":
                    fileUpdate.write(f"{ID}\t{covarDict[ID][covar]}\n")
                fileOut.write(f"\t{covarDict[ID][covar]}")
            for PCNum in range(1, i):
                PC = f"PC{PCNum}"
                fileOut.write(f"\t{dictPCA[ID][PC]}")
            fileOut.write("\n")
    fileOut.close()
    fileUpdate.close()

    if "stepwise" in config['model']:
        R = config['programs']['r']
        stepwise = config['model']['stepwise']
        pheno = config["input"]["pheno"].upper()

        data = pd.read_table(f"{folder}/{name}.tsv")
        data[pheno] = data[pheno] -1
        data.to_csv(f"{folder}/{name}_Step.tsv", index = False, sep = "\t")

        execute(f"{R} {stepwise} {folder}/{name}_Step.tsv {pheno} {folder}/{name}", logFile, debug, True)
        fileIn = open(f"{folder}/{name}_variables.tsv")

        covarReturn = ""
        for line in fileIn:
            if covarReturn == "":
                covarReturn = line.strip()
            else:
                covarReturn = f"{covarReturn} {line.strip()}"

        return f"{folder}/{name}.tsv", pheno, covarReturn, f"{folder}/{name}_UpdateSex"

    elif "finalmodel" in config['model']:
        pheno, covar = config['model'].split("~")
        return f"{folder}/{name}.tsv", pheno, covar.replace("+", " "), f"{folder}/{name}_UpdateSex"
    else:
        sys.exit("Error: Model has no valid option (script or regression)")
