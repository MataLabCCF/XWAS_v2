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


def getPCANonProjected(bfile, runPCA, folder, plink1, logFile):
    print(f"We will run X-PCA without reference")
    PCA = getPCA(bfile, runPCA, "PCA_Autosomal", folder, plink1, logFile)
    return PCA


def getXPCANonProjected(bfileFemale, bfileMale, runPCA, folder, plink1, logFile):
    print(f"We will run X-PCA without reference")
    XPCAFemale = getPCA(bfileFemale, runPCA, "XPCA_Female", folder, plink1, logFile)
    XPCAMale = getPCA(bfileMale, runPCA, "XPCA_Male", folder, plink1, logFile)

    return XPCAFemale, XPCAMale


def mergeAndGetProjectedPCA(bfileTarget, bfileRef, folder, name, type, plink1, gcta, threads, X, logFile):
    nameOut = f"{name}_{type}"
    bfileMerged, targetCommon, refCommon = mergeRefAndTarget(bfileTarget, bfileRef, folder, nameOut, plink1, logFile)
    targetLD, refLD = removeLDAndMAFToPCA(bfileMerged, targetCommon, refCommon, folder, nameOut, plink1, logFile)
    return getProjectedPCA(targetLD, refLD, gcta, X, threads, type, folder, plink1, logFile)


def getXPCAProjected(XRef, bfileFemale, bfileMale, gcta, folder, name, plink1, threads, logFile):
    print(f"Spliting the reference data {XRef} using the PLINK impute-sex")
    bfileMaleRef, bfileFemaleRef = separateGenotypedDataBySexWithoutCovar(XRef, f"{name}_RefX", folder, plink1, logFile)
    XPCAFemale = mergeAndGetProjectedPCA(bfileFemale, bfileFemaleRef, folder, name, "XPCA_Female", plink1, gcta,
                                         threads, True, logFile)
    XPCAMale = mergeAndGetProjectedPCA(bfileMale, bfileMaleRef, folder, name, "XPCA_Male", plink1, gcta, threads, True,
                                       logFile)

    return XPCAFemale, XPCAMale


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
                      action="store_true", default=False)
    data.add_argument('--threads', help='Number of processors to be used (default = 1)',
                      required=False, default=1)

    refData = parser.add_argument_group("Reference data arguments")
    refData.add_argument('-r', '--XRef', help='Data with parental reference data to run XPCA (optional)',
                         required=False, default="")
    refData.add_argument('-R', '--AutosomalRef',
                         help='Data with parental reference data to run autosomal PCA (optional)', required=False,
                         default="")

    output = parser.add_argument_group("Output arguments")
    output.add_argument('-n', '--name', help='Analysis name', required=True)
    output.add_argument('-f', '--folder', help='Folder to output files', required=True)

    pcaArg = parser.add_argument_group("PCA arguments")
    pcaArg.add_argument('--model', help='Regression model. If you do not provide a model the script will'
                                        'use the selectModel.R to build the model', required=False, default="",
                        nargs="+")
    pcaArg.add_argument('--PCA', help='PCA to be calculated and used. (default = 1) \n'
                                      '1- X-PCA\n2- Projected X-PCA\n3- Autosomal PCA\n4- Projected autosomal PCA\n'
                                      '5- X-PCA + Autosomal PA\n 6- Projected X-PAS + Projected Autosomal\n',
                        required=False, default="1")

    regressionArg = parser.add_argument_group("Regression arguments")
    regressionArg.add_argument('--rsquare', help='R2 imputation cutoff to regression (default = 0.8)', required=False,
                               default=0.8)
    regressionArg.add_argument('--firth', help='Force PLINK2 firth regression (default Plink2 will decide)',
                               required=False, default=False, action="store_true")

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

    if args.PCA == "1" or args.PCA == "5":
        XPCAFemale, XPCAMale = getXPCANonProjected(bfileFemale, bfileMale, args.runPCA, args.folder, args.plink1,
                                                   logFile)
        PCASource.append("XPCA")

        covarDict, PCList = addPCAToCovarDict(XPCAFemale, covarDict, "XPCA")
        covarDict, PCList = addPCAToCovarDict(XPCAMale, covarDict, "XPCA")
    if args.PCA == "2" or args.PCA == "6":
        if args.XRef != "":
            XPCAFemale, XPCAMale = getXPCAProjected(args.XRef, bfileFemale, bfileMale, args.gcta, args.folder,
                                                    args.name, args.plink1, args.threads, logFile)
            PCASource.append("XPCA_Projected")

            covarDict, PCList = addPCAToCovarDict(XPCAFemale, covarDict, "XPCA_Projected")
            covarDict, PCList = addPCAToCovarDict(XPCAMale, covarDict, "XPCA_Projected")
        else:
            finish(f"The PCA option 2 requires the XRef file")

    if args.PCA == "3" or args.PCA == "5":
        autosomalPCA = getPCANonProjected(args.autosomal, args.runPCA, args.folder, args.plink1, logFile)

        PCASource.append("AutosomalPCA")
        covarDict, PCList = addPCAToCovarDict(autosomalPCA, covarDict, "AutosomalPCA")

    if args.PCA == "4" or args.PCA == "6":
        if args.AutosomalRef != "":
            autosomalPCA = mergeAndGetProjectedPCA(args.autosomal, args.AutosomalRef, args.folder, args.name,
                                                   "PCA_Projected", args.plink1, args.gcta, args.threads,
                                                   False, logFile)

            PCASource.append("AutosomalPCA_Projected")
            covarDict, PCList = addPCAToCovarDict(autosomalPCA, covarDict, "AutosomalPCA_Projected")
        else:
            finish(f"The PCA option 4 requires the AutosomalRef file")


    if args.model:
        PCList = []
        for var in args.model:
            modelMale.append(var)
            modelFemale.append(var)
            if "PC" in var:
                PCList.append(var)

    femaleOutlierList = createOutlierList(covarDict, bfileFemale, PCASource, PCList, args.folder, f"{args.name}_Female")
    maleOutlierList = createOutlierList(covarDict, bfileMale, PCASource, PCList, args.folder, f"{args.name}_Male")

    covarMale = buildCovarFile(covarDict, maleOutlierList, PCList, PCASource, bfileMale, args.folder,
                               f"{args.name}_covarMale")
    covarFemale = buildCovarFile(covarDict, femaleOutlierList, PCList, PCASource, bfileFemale, args.folder,
                                 f"{args.name}_covarFemale")

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
