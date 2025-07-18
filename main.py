import os
import sys
import argparse
from Analysis.PCA import *
from Analysis.Regression import *
from Handlers.COVAR import *
from Handlers.PLINK1 import *
from Handlers.PLINK2 import *

def detectOutliersBoth(PCAModel, bfileX, bfileAuto, XRef, autosomalRef, gcta, plink1, folder, threads, logFile):
    outlierList = []

    print("===========================================================================================")
    print(f"Removing outliers for both to PCAModel {PCAModel} on detectOutliersBoth")


    if PCAModel == "1" or PCAModel == "5":
        print(f"\tRemoving LD for female non-project XPCA")
        femaleNoLD = removeLD(bfileX, "NonProjected_XPCA", folder, plink1, logFile)
        XPCA = getPCA(femaleNoLD, "", gcta, True, threads, "NonProjected_XPCA_Both", folder, plink1, logFile)
        outlierList = getOutlierFromPCA(XPCA, outlierList)

    # Getting outlier from projected X PCA
    if PCAModel == "2" or PCAModel == "6":
        if XRef != "":
            XPCA = mergeAndGetProjectedPCA(bfileX, XRef, folder, "Projected", "XPCA_Both", plink1, gcta, threads, True, logFile)
            outlierList = getOutlierFromPCA(XPCA, outlierList)
        else:
            finish(f"The PCA option {PCAModel} requires the XRef file", logFile)

    if PCAModel == "3" or PCAModel == "5":
        autosomalNoLD = removeLD(bfileAuto, "NonProjected_AutosomalPCA", folder, plink1, logFile)
        PCA = getPCA(autosomalNoLD, "", gcta, False, threads, "NonProjected_Auto_Both", folder, plink1, logFile)
        outlierList = getOutlierFromPCA(PCA, outlierList)

    if PCAModel == "4" or PCAModel == "6":
        if args.AutosomalRef != "":
            PCA = mergeAndGetProjectedPCA(bfileAuto, autosomalRef, folder, "Projected","AutosomalPCA", plink1,
                                          gcta, threads,False, logFile)
            outlierList = getOutlierFromPCA(PCA, outlierList)

        else:
            finish(f"The PCA option {PCAModel} requires the AutosomalRef file", logFile)

    return outlierList

def detectOutliers(PCAModel, bfileFemale, bfileMale, autosomalFemale, autosomalMale, XRef, autosomalRef, gcta, plink1, folder, threads, logFile):
    outlierList = []
    print(f"Removing outliers")

    bfileFemaleRef = bfileMaleRef = ""

    #Getting outlier from non projected X PCA
    if PCAModel == "1" or PCAModel == "5":
        print(f"\tRemoving LD for female non-project XPCA")
        femaleNoLD = removeLD(bfileFemale, "NonProjected_XPCA", folder, plink1, logFile)
        XPCA = getPCA(femaleNoLD, "", gcta, True, threads, "NonProjected_XPCA_Female", folder, plink1, logFile)
        outlierList = getOutlierFromPCA(XPCA, outlierList)

        maleNoLD = removeLD(bfileMale, "NonProjected_XPCA", folder, plink1, logFile)
        XPCA = getPCA(maleNoLD, "", gcta, True, threads, "NonProjected_XPCA_Male", folder, plink1, logFile)
        outlierList = getOutlierFromPCA(XPCA, outlierList)


    # Getting outlier from projected X PCA
    if PCAModel == "2" or PCAModel == "6":
        if XRef != "":
            print(f"Spliting the reference data {XRef} using the PLINK impute-sex")
            bfileMaleRef, bfileFemaleRef = separateGenotypedDataBySexWithoutCovar(XRef, f"RefX", folder, plink1, logFile)
            XPCA = mergeAndGetProjectedPCA(bfileFemale, bfileFemaleRef, folder, "Projected", "XPCA_Female", plink1, gcta,
                                                 threads, True, logFile)
            outlierList = getOutlierFromPCA(XPCA, outlierList)


            XPCA = mergeAndGetProjectedPCA(bfileMale, bfileMaleRef, folder, "Projected", "XPCA_Male", plink1, gcta,
                                               threads, True, logFile)
            outlierList = getOutlierFromPCA(XPCA, outlierList)

        else:
            finish(f"The PCA option 2 and 6 requires the XRef file", logfile)

    # Getting outlier from non projected autosomal PCA
    if PCAModel == "3" or PCAModel == "5":
        autosomalNoLD = removeLD(autosomalFemale, "NonProjected_AutosomalPCA", folder, plink1, logFile)
        PCA = getPCA(autosomalNoLD, "", gcta, False , threads, "NonProjected_Auto", folder, plink1, logFile)
        outlierList = getOutlierFromPCA(PCA, outlierList)

        autosomalNoLD = removeLD(autosomalMale, "NonProjected_AutosomalPCA", folder, plink1, logFile)
        PCA = getPCA(autosomalNoLD, "", gcta, False, threads, "NonProjected_Auto", folder, plink1, logFile)
        outlierList = getOutlierFromPCA(PCA, outlierList)

    # Getting outlier from projected autosomal PCA
    if PCAModel == "4" or PCAModel == "6":
        if args.AutosomalRef != "":
            PCA = mergeAndGetProjectedPCA(autosomalMale, autosomalRef, folder, "Projected","AutosomalPCA", plink1,
                                          gcta, threads,False, logFile)
            outlierList = getOutlierFromPCA(PCA, outlierList)

            PCA = mergeAndGetProjectedPCA(autosomalFemale, autosomalRef, folder, "Projected", "AutosomalPCA", plink1,
                                          gcta, threads, False, logFile)
            outlierList = getOutlierFromPCA(PCA, outlierList)

        else:
            finish(f"The PCA option 4 requires the AutosomalRef file", logfile)
    return outlierList

def createFolder(folder, log=True):
    if not os.path.exists(folder):
        os.mkdir(folder)
    if log:
        logFile = open(f"{folder}/XWAS.log", "w")
        return logFile

def finish(message, logFile):
    logFile.write(message)
    sys.exit(message)

def removeOutlierAndMissing(bfile, covarDict, outlierList, nonOutlierFolder, outName, plink1, logFile):
    fam = open(f"{bfile}.fam")
    samplesSelected = open(f"{nonOutlierFolder}/{outName}_toKeep.txt", "w")
    header = True

    sampleList = []
    for line in fam:
        if header:
            header = False
        else:
            IID = line.strip().split()[1]
            if IID not in outlierList:
                remove = False
                for covar in covarDict[IID]:
                    if covarDict[IID][covar] == "" or covarDict[IID][covar] == "NA" or covarDict[IID][covar] == "-9":
                        print(f"{IID} -> {covar} --> {covarDict[IID][covar]}")
                        remove = True

                if not remove:
                    sampleList.append(IID)
                    samplesSelected.write(f"{IID}\t{IID}\n")
    samplesSelected.close()
    fam.close()
    return extractSamples(bfile, f"{nonOutlierFolder}/{outName}_toKeep.txt", outName,
                                        nonOutlierFolder, plink1, logFile)


def mergeAndGetProjectedPCA(bfileTarget, bfileRef, folder, name, type, plink1, gcta, threads, X, logFile):
    nameOut = f"{name}_{type}"


    bfileMerged, targetCommon, refCommon = mergeRefAndTarget(bfileTarget, bfileRef, folder, nameOut, plink1, logFile)
    targetLD, refLD = removeLDAndMAFToPCA(bfileMerged, targetCommon, refCommon, folder, nameOut, plink1, logFile)
    return getPCA(targetLD, refLD, gcta, X, threads, type, folder, plink1, logFile)



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
    data.add_argument('-c', '--continuous', help='Set as continuous phenotype (use beta, not OR)', required=False,
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
    pcaArg.add_argument('--model', help='Regression model. List of covariates to be used separated by comma '
                                        '(do not add spaces). If not used, we will use stepwise AIC', required=False, default="")
    pcaArg.add_argument('--modelSelection', help=' Path to modelSelection.R script provided on the XWAS_v2', required=False,
                        default="")
    pcaArg.add_argument('--PCAOutlier', help='PCA to be calculated to outlier removal. (default = 1) \n'
                                      '1- X-PCA\n2- Projected X-PCA\n3- Autosomal PCA\n4- Projected autosomal PCA\n'
                                      '5- X-PCA + Autosomal PCA\n 6- Projected X-PCA + Projected Autosomal\n',
                        required=False, default="1")
    pcaArg.add_argument('--PCARegression', help='PCA to be calculated to regression. (default = 1) \n'
                                             '1- X-PCA\n2- Autosomal PCA\n',
                        required=False, default="1")

    regressionArg = parser.add_argument_group("Regression arguments")
    regressionArg.add_argument('--rsquare', help='R2 imputation cutoff to regression (default = 0.8)', required=False,
                               default=0.8)
    regressionArg.add_argument('--firth', help='Force PLINK2 firth regression (default Plink2 will decide)',
                               required=False, default=False, action="store_true")
    regressionArg.add_argument('--phenoName', help='Name of the phenotype to be tested on covar file (default: DISEASE)',
                               required=False, default="DISEASE")

    programs = parser.add_argument_group("Programs")
    programs.add_argument('--gwama', help='GWAMA program (default = gwama)', required=False, default="gwama")
    programs.add_argument('--plink2', help='Path of PLINK 2 (default = plink2)', required=False, default="plink2")
    programs.add_argument('--plink1', help='Path of PLINK 1 (default = plink)', required=False, default="plink")
    #programs.add_argument('--runPCA', help='Path of runPCA script', required=True)
    programs.add_argument('--gcta', help='Path of gcta', required=True)

    args = parser.parse_args()

    logFile = createFolder(args.folder)

    debug = False
    if not debug:
        covarDict = readCovarFile(args.tableCovar)

        print(f"Spliting the target data {args.chrX} using the covar table")
        bfileFemale, bfileMale, bfileBoth = separateGenotypedDataBySex(args.chrX, covarDict, f"{args.name}_TargetX",
                                                            args.folder, args.plink1, logFile)

        bfileAutoFemale = bfileAutoMale = bfileAutoBoth = ""
        if args.autosomal:
            bfileAutoFemale, bfileAutoMale, bfileAutoBoth = separateGenotypedDataBySex(args.autosomal, covarDict, f"{args.name}_TargetAutosomal",
                                                                        args.folder, args.plink1, logFile)

        #Calculate PCA to remove outliers

        outlierFolder = f"{args.folder}/Outliers/"
        createFolder(f"{outlierFolder}")
        outlierList = detectOutliers(args.PCAOutlier, bfileFemale, bfileMale, bfileAutoFemale, bfileAutoMale,
                                     args.XRef, args.AutosomalRef, args.gcta, args.plink1, outlierFolder, args.threads, logFile)



        # We will not use projected PCA to GWAS analysis because we do not have good references for our population. Here we will implement X-PCA or autosomal PCA
        nonOutlierFolder = f"{args.folder}/NonOutliers/"
        createFolder(f"{nonOutlierFolder}")
        folderPCA = f"{args.folder}/PCARegression/"
        createFolder(f"{folderPCA}")
        if args.PCARegression == "1":
            # My outliers are now a list. Let's remove the outliers, generate the PCs to regression, perform the regression
            # and meta-analysis
            bfileFinalFemale = removeOutlierAndMissing(bfileFemale, covarDict, outlierList, nonOutlierFolder,
                                                       "Female", args.plink1, logFile)
            bfileFinalMale = removeOutlierAndMissing(bfileMale ,covarDict, outlierList, nonOutlierFolder, "Male",
                                                     args.plink1, logFile)

            femaleNoLD = removeLD(bfileFinalFemale, "NonProjected_XPCA", folderPCA, args.plink1, logFile)
            XPCAFemale = getPCA(femaleNoLD, "", args.gcta, True, args.threads, "NonProjected_XPCA_Female",
                          folderPCA, args.plink1, logFile)
            maleNoLD = removeLD(bfileFinalMale, "NonProjected_XPCA", folderPCA, args.plink1, logFile)
            XPCAMale = getPCA(maleNoLD, "", args.gcta, True, args.threads, "NonProjected_XPCA_Male",
                          folderPCA, args.plink1, logFile)

            covarMale = createCovarFileToRegression(covarDict, XPCAMale, folderPCA, "CovarMale")
            covarFemale = createCovarFileToRegression(covarDict, XPCAFemale, folderPCA, "CovarFemale")

        elif args.PCARegression == "2":
            bfileFinalFemale = removeOutlierAndMissing(bfileAutoFemale, covarDict, outlierList, nonOutlierFolder,
                                                       "Female", args.plink1, logFile)
            bfileFinalMale = removeOutlierAndMissing(bfileAutoMale, covarDict, outlierList, nonOutlierFolder, "Male",
                                                     args.plink1, logFile)

            femaleNoLD = removeLD(bfileFinalFemale, "NonProjected_AutoPCA", folderPCA, args.plink1, logFile)
            PCAFemale = getPCA(femaleNoLD, "", args.gcta, True, args.threads, "NonProjected_AutoPCA_Female",
                                folderPCA, args.plink1, logFile)
            maleNoLD = removeLD(bfileFinalMale, "NonProjected_AutoPCA", folderPCA, args.plink1, logFile)
            PCAMale = getPCA(maleNoLD, "", args.gcta, True, args.threads, "NonProjected_AutoPCA_Male",
                              folderPCA, args.plink1, logFile)

            covarMale = createCovarFileToRegression(covarDict, PCAMale, folderPCA, "CovarMale")
            covarFemale = createCovarFileToRegression(covarDict, PCAFemale, folderPCA, "CovarFemale")

        else:
            exit(f"Invalid option to --PCARegression")

        folderReg = f"{args.folder}/Regression/"
        createFolder(f"{folderReg}")
        createFolder(f"{folderReg}/Results")

        pvarFileMale = convertVCFToPFile(args.imputed, bfileFinalMale, covarMale,f"{args.name}_toRegressionMale",
                                         folderReg, args.plink2, logFile)
        pvarFileFemale = convertVCFToPFile(args.imputed, bfileFinalFemale, covarFemale, f"{args.name}_toRegressionFemale",
                                           folderReg, args.plink2, logFile)


        #=================================================== FOR BOTH ===================================================
        #Generating the both dataset
        bfileBothNonOutlier =  removeOutlierAndMissing(bfileBoth, covarDict, outlierList, nonOutlierFolder,
                                                       "BothX", args.plink1, logFile)
        bfileAutoBothNonOutlier = ""
        if args.autosomal:
            bfileAutoBothNonOutlier = removeOutlierAndMissing(bfileAutoBoth, covarDict, outlierList, nonOutlierFolder,
                                                              "BothAuto", args.plink1, logFile)



        outlierListBoth = detectOutliersBoth(args.PCAOutlier, bfileBothNonOutlier, bfileAutoBothNonOutlier, args.XRef,
                                             args.AutosomalRef, args.gcta, args.plink1, outlierFolder, args.threads, logFile)

        if args.PCARegression == "1":
            bfileFinalBoth = removeOutlierAndMissing(bfileBothNonOutlier, covarDict, outlierListBoth, nonOutlierFolder,
                                                       "Both", args.plink1, logFile)
            bothNoLD = removeLD(bfileFinalBoth, "NonProjected_XPCA", folderPCA, args.plink1, logFile)
            XPCABoth = getPCA(bothNoLD, "", args.gcta, True, args.threads, "NonProjected_XPCA_Both",
                                folderPCA, args.plink1, logFile)
            covarBoth = createCovarFileToRegression(covarDict, XPCABoth, folderPCA, "CovarBoth")
        elif args.PCARegression == "2":
            bfileFinalBoth = removeOutlierAndMissing(bfileAutoBothNonOutlier, covarDict, outlierList, nonOutlierFolder,
                                                       "Both", args.plink1, logFile)
            bothNoLD = removeLD(bfileFinalBoth, "NonProjected_AutoPCA", folderPCA, args.plink1, logFile)
            PCABoth = getPCA(bothNoLD, "", args.gcta, False, args.threads, "NonProjected_AutoPCA_Both",
                               folderPCA, args.plink1, logFile)
            covarBoth = createCovarFileToRegression(covarDict, PCABoth, folderPCA, "CovarBoth")

        pvarFileBoth = convertVCFToPFile(args.imputed, bfileFinalBoth, covarBoth, f"{args.name}_toRegressionBoth",
                                         folderReg, args.plink2, logFile)


    print("I escaped all preparation steps")
    folderReg = f"{args.folder}/Regression/"
    pvarFileMale = f"{folderReg}{args.name}_toRegressionMale"
    pvarFileFemale = f"{folderReg}{args.name}_toRegressionFemale"

    folderPCA = f"{args.folder}/PCARegression/"
    covarMale = f"{folderPCA}CovarMale.tsv"
    covarFemale = f"{folderPCA}CovarFemale.tsv"


    regressionMale = runRegressionPlink2(pvarFileMale, covarMale, args.model, args.rsquare, args.firth, args.phenoName,
                                         args.threads, args.plink2, folderReg, f"{args.name}_RegMale", args.modelSelection, logFile)
    regressionFemale = runRegressionPlink2(pvarFileFemale, covarFemale, args.model, args.rsquare, args.firth, args.phenoName,
                                           args.threads, args.plink2, folderReg, f"{args.name}_RegFemale", args.modelSelection,logFile)
    regressionBoth = runRegressionPlink2(pvarFileBoth, covarBoth, args.model, args.rsquare, args.firth, args.phenoName,
                                           args.threads, args.plink2, folderReg, f"{args.name}_RegBoth", args.modelSelection,logFile)

    gwamaMetaAnalysis(regressionFemale, regressionMale, pvarFileFemale, pvarFileMale, args.firth,
                      args.plink2, args.name, folderReg, args.gwama, logFile)






