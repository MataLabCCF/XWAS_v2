import os
import sys
import argparse
from Analysis.PCA import *
from Analysis.Regression import *
from Handlers.COVAR import *
from Handlers.PLINK1 import *
from Handlers.PLINK2 import *

def createFolder(folder, log=True):
    if not os.path.exists(folder):
        os.mkdir(folder)
    if log:
        logFile = open(f"{folder}/XWAS.log", "w")
        return logFile


def finish(message, logFile):
    logFile.write(message)
    sys.exit(message)



# Differences from XWAS v1 and XWAS v2
# - We removed the country file. It was implemented on v1 and never used because the analysis per country
# requires a QC per country.


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PCA and regression')

    data = parser.add_argument_group("Data arguments")
    data.add_argument('-X', '--chrX', help='Genotyped file name for X chromosome', required=False)
    data.add_argument('-A', '--autosomal', help='Genotyped file name for autosomal chromosome', required=False)
    data.add_argument('-t', '--tableCovar', help='File with covariatives to be added to the model', required=True)
    data.add_argument('-i', '--imputed', help='Imputed file name', required=False)
    data.add_argument('-c', '--continuos', help='Set as continuos phenotype (use beta, not OR)', required=False,
                      action = "store_true", default = False)
    data.add_argument('--threads', help='Number of processors to be used (default = 1)',
                      required=False, default = 1)

    refData = parser.add_argument_group("Reference data arguments")
    refData.add_argument('-r', '--XRef', help='Data with parental reference data to run XPCA (optional)', required=False, default="")
    refData.add_argument('-R', '--AutosomalRef', help='Data with parental reference data to run autosomal PCA (optional)', required=False, default="")

    output = parser.add_argument_group("Output arguments")
    output.add_argument('-n', '--name', help='Analysis name', required=True)
    output.add_argument('-f', '--folder', help='Folder to output files', required=True)

    pcaArg = parser.add_argument_group("PCA arguments")
    pcaArg.add_argument('--model', help='Regression model. If you do not provide a model the script will'
                                        'use the selectModel.R to build the model', required=False, default="", nargs="+")

    regressionArg = parser.add_argument_group("Regression arguments")
    regressionArg.add_argument('--rsquare', help='R2 imputation cutoff to regression (default = 0.8)', required=False, default= 0.8)
    regressionArg.add_argument('--firth', help='Force PLINK2 firth regression (default Plink2 will decide)',
                               required=False, default = False, action = "store_true")

    programs = parser.add_argument_group("Programs")
    programs.add_argument('--gwama', help='GWAMA program (default = gwama)', required=False, default="gwama")
    programs.add_argument('--plink2', help='Path of PLINK 2 (default = plink2)', required=False, default="plink2")
    programs.add_argument('--plink1', help='Path of PLINK 1 (default = plink)', required=False, default="plink")
    programs.add_argument('--runPCA', help='Path of runPCA script', required=True)
    programs.add_argument('--gcta', help='Path of gcta', required=True)
    programs.add_argument('--selectModel', help='Path of selectModel script', required=False)

    args = parser.parse_args()

    logFile = createFolder(args.folder)

    covarDict = readCovarFile(args.tableCovar)

    modelMale = []
    modelFemale = []

    print(f"Spliting the target data {args.chrX} using the covar table")
    bfileFemale, bfileMale = separateGenotypedDataBySex(args.chrX, covarDict, f"{args.name}_TargetX",
                                                        args.folder, args.plink1, logFile)

    PCASource = []

    if args.XRef == "" and args.AutosomalRef == "":
        print(f"We will run X-PCA without reference")
        XPCAFemale = getPCA(bfileFemale, args.runPCA, "XPCA_Female", args.folder, args.plink1, logFile)
        XPCAMale = getPCA(bfileMale, args.runPCA, "XPCA_Male", args.folder, args.plink1, logFile)

        PCASource.append("XPCA_NR")

        covarDict, PCList = addPCAToCovarDict(XPCAFemale, covarDict, "XPCA_NR")
        covarDict, PCList = addPCAToCovarDict(XPCAMale, covarDict, "XPCA_NR")


    if args.XRef != "":
        print(f"Spliting the reference data {args.XRef} using the PLINK impute-sex")
        bfileMaleRef, bfileFemaleRef = separateGenotypedDataBySexWithoutCovar(args.XRef, f"{args.name}_RefX",
                                                                              args.folder, args.plink1, logFile)

        bfileFemaleMerged, targetFemaleCommon, refFemaleCommon = mergeRefAndTarget(bfileFemale, bfileFemaleRef, args.folder,f"{args.name}_Female_XPCA", args.plink1, logFile)
        targetFemaleLD, refFemaleLD = removeLDAndMAFToPCA(bfileFemaleMerged, targetFemaleCommon, refFemaleCommon, args.folder,f"{args.name}_Female_XPCA", args.plink1, logFile)

        bfileMaleMerged, targetMaleCommon, refMaleCommon = mergeRefAndTarget(bfileMale, bfileMaleRef, args.folder,f"{args.name}_Male_LD", args.plink1, logFile)
        targetMaleLD, refMaleLD = removeLDAndMAFToPCA(bfileMaleMerged, targetMaleCommon, refMaleCommon, args.folder, f"{args.name}_Male_LD", args.plink1, logFile)

        XPCAFemale = getProjectedPCA(targetFemaleCommon, refFemaleCommon, args.gcta, True, args.threads, "XPCA_Female", args.folder, args.plink1, logFile)
        XPCAMale = getProjectedPCA(targetMaleCommon, refMaleCommon, args.gcta, True, args.threads, "XPCA_Male", args.folder, args.plink1, logFile)

        PCASource.append("XPCA_WR")

        covarDict, PCList = addPCAToCovarDict(XPCAFemale, covarDict, "XPCA_WR")
        covarDict, PCList = addPCAToCovarDict(XPCAMale, covarDict, "XPCA_WR")

    if args.AutosomalRef != "":
        bfileMerged, targetCommon, refCommon = mergeRefAndTarget(args.autosomal, args.AutosomalRef, args.folder, f"{args.name}_AutosomalPCA", args.plink1, logFile)
        targetLD, refLD = removeLDAndMAFToPCA(bfileMerged, targetCommon, refCommon, args.folder, f"{args.name}_Male_LD", args.plink1, logFile)

        autosomalPCA = getProjectedPCA(targetCommon, refCommon, args.gcta, False, args.threads, "AutosomalPCA", args.folder, args.plink1, logFile)

        PCASource.append("Autosomal_WR")

        covarDict, PCList = addPCAToCovarDict(autosomalPCA, covarDict, "Autosomal_WR")

    if args.model:
        PCList = []
        for var in args.model:
            modelMale.append(var)
            modelFemale.append(var)
            if "PC" in var:
                PCList.append(var)

    femaleOutlierList = createOutlierList(covarDict, bfileFemale, PCASource, PCList, args.folder, f"{args.name}_Female")
    maleOutlierList = createOutlierList(covarDict, bfileMale, PCASource, PCList, args.folder, f"{args.name}_Male")

    covarMale = buildCovarFile(covarDict, maleOutlierList, PCList, PCASource, bfileMale, args.folder, f"{args.name}_covarMale")
    covarFemale = buildCovarFile(covarDict, femaleOutlierList, PCList, PCASource, bfileFemale, args.folder, f"{args.name}_covarFemale")

    if not args.model:
        modelMale = modelSelection(covarMale, args.selectModel, args.folder, f"{args.name}_Male_model", logFile)
        modelFemale = modelSelection(covarFemale, args.selectModel, args.folder, f"{args.name}_Female_model", logFile)

    pvarFileMale = convertVCFToPFile(args.imputed, bfileMale, maleOutlierList, covarDict,
                                     f"{args.name}_toRegressionMale", args.folder, args.plink2, logFile)
    pvarFileFemale = convertVCFToPFile(args.imputed, bfileFemale, femaleOutlierList, covarDict,
                                       f"{args.name}_toRegressionFemale", args.folder, args.plink2, logFile)

    regressionMale = runRegressionPlink2(pvarFileMale, covarMale, modelMale, args.rsquare, args.firth, args.threads,
                                         args.plink2, args.folder, f"{args.name}_RegMale", logFile)
    regressionFemale = runRegressionPlink2(pvarFileFemale, covarFemale, modelFemale, args.rsquare, args.firth,
                                           args.threads, args.plink2, args.folder, f"{args.name}_RegFemale", logFile)


    gwamaMetaAnalysis(regressionFemale, regressionMale, pvarFileFemale, pvarFileMale, args.firth,
                                       args.plink2, args.name, args.folder, args.gwama, logFile)