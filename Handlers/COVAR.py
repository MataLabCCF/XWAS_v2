


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