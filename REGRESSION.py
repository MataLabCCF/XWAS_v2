from UTILS import *


def prepareInputGWAMA(regression, pfile, isFirth, plink2, name, folder):
    dictID = {}
    pvarFile = open(f'{pfile}.pvar')

    #Info from VCF-like PVAR file
    header = True
    for line in pvarFile:
        if header:
            if "#CHROM" in line:
                header = False
        else:
            split = line.strip().split()

            CHROM = split[0]
            POS = split[1]
            ID = split[2]
            REF = split[3]
            ALT = split[4]

            dictID[ID] = {}
            dictID[ID]["newID"] = f"{CHROM}:{POS}:{REF}:{ALT}"
            dictID[ID]["imputed"] = 1

            if "TYPED" in line:
                dictID[ID]["imputed"] = 0
    pvarFile.close()

    os.system(f"{plink2} --pfile {pfile} --freq --out {folder}/{name}_freq")
    file = open(f"{folder}/{name}_freq.afreq")

    header = True
    for line in file:
        if header:
            header = False
        else:
            CHROM, ID, REF, ALT, ALT_FREQS, OBS_CT = line.strip().split()
            dictID[ID]["ALT"] = ALT
            dictID[ID]["ALT_FREQ"] = ALT_FREQS
            dictID[ID]["N"] = OBS_CT
    file.close()


    regressionName = f"{regression}.DISEASE.glm.logistic.hybrid"
    if isFirth:
        regressionName = f"{regression}.DISEASE.glm.firth"

    file = open(regressionName)
    fileGWAMA = open(f"{folder}/{name}.in", "w")
    fileGWAMA.write("MARKERNAME\tCHR\tPOS\tIMPUTED\tN\tEA\tNEA\tEAF\tOR\tOR_95L\tOR_95U\n")

    header = True
    for line in file:
        if header:
            dictColHeader = {}

            split = line.strip().split()
            interest = ["#CHROM", "POS", "ID", "REF", "ALT", "A1", "OR", "LOG(OR)_SE", "P", "OBS_CT", "U95", "L95"]

            for i in range(0, len(split)):
                if split[i] in interest:
                    #print(f"Col {split[i]} -> index {i}")
                    dictColHeader[split[i]] = i
            header = False
        else:
            split = line.strip().split()
            if split[dictColHeader["OR"]] != "NA":
                ID = split[dictColHeader["ID"]]
                CHR = split[dictColHeader["#CHROM"]]
                POS = split[dictColHeader["POS"]]
                IMPUTED = dictID[ID]["imputed"]
                N = split[dictColHeader["OBS_CT"]]
                EA = split[dictColHeader["A1"]]
                if EA == split[dictColHeader["ALT"]]:
                    NEA = split[dictColHeader["REF"]]
                else:
                    NEA = split[dictColHeader["ALT"]]

                if EA == dictID[ID]["ALT"]:
                    EAF = dictID[ID]["ALT_FREQ"]
                else:
                    EAF = 1 - float(dictID[ID]["ALT_FREQ"])

                OR = split[dictColHeader["OR"]]
                OR_95L = split[dictColHeader["L95"]]
                OR_95U = split[dictColHeader["U95"]]

                fileGWAMA.write(f"{ID}\t{CHR}\t{POS}\t{IMPUTED}\t{N}\t{EA}\t{NEA}\t{EAF}\t{OR}\t{OR_95L}\t{OR_95U}\n")
    fileGWAMA.close()
    return f"{folder}/{name}.in"

def gwamaMetaAnalysis(regressionFemale, regressionMale, pvarFileFemale, pvarFileMale, isFirth, name, folder, config, logFile, debug):
    plink2 = config["programs"]["plink2"]
    gwama = config["programs"]["gwama"]

    inputFemale = prepareInputGWAMA(regressionFemale, pvarFileFemale, isFirth, plink2, f"{name}_FemaleGWAMA", folder)
    inputMale = prepareInputGWAMA(regressionMale, pvarFileMale, isFirth, plink2, f"{name}_MaleGWAMA", folder)

    inputGWAMA = open(f"{folder}/inputGWAMA.in", "w")
    inputGWAMA.write(f"{inputFemale}\tF\n")
    inputGWAMA.write(f"{inputMale}\tM\n")
    inputGWAMA.close()

    execute(f"{gwama} -i {folder}/inputGWAMA.in -o {folder}/{name} -r -gc --sex", logFile, debug)

def runRegressionPlink2(pvar, covarFile, rsquare, firth, phenoName, covars, folder, name, config, threads, logFile, debug):
    plink2 = config["programs"]["plink2"]


    outputFolder = f"{folder}/Results/"
    createFolder(f"{outputFolder}", False)

    outputPrefix = f"{outputFolder}/{name}"

    command = (f"{plink2} --pfile {pvar}  --covar-variance-standardize --ci 0.95 "
               f"--covar {covarFile} --covar-name {covars} --out {outputPrefix} "
               f"--pheno {covarFile} --pheno-name {phenoName} --threads {threads} --glm hide-covar")

    if firth:
        command = command + " firth"
    if rsquare != "":
        command = command + f" --extract-if-info \"R2 > {rsquare}\""
    execute(command, logFile, debug)

    return outputPrefix