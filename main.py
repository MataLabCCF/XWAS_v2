import argparse
from PCA import *
from UTILS import *
from PLINK import *
from REGRESSION import *


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mata Lab XWAS v2.1')

    data = parser.add_argument_group("Input arguments")
    data.add_argument('-c', '--config', help='Path to configuration file', required=True)
    data.add_argument('-t', '--threads', help='Number of threads to be used by PLINK and GCTA (default = 1)', required=False, default= 1)
    data.add_argument('-m', '--mem', help='Memory parameter used by PLINK (default = None, Plink will use 50% of computer RAM memory', required=False, default="")

    out = parser.add_argument_group("Output arguments")
    out.add_argument('-f', '--folder', help='Folder name', required=True)

    args = parser.parse_args()

    print(f"Reading the configuration file")

    config, covarDict = readConfigFile(args.config)

    logFile = createFolder(args.folder)

    debug = False

    #Step 1: Split the typed data by males and females
    print(f"Spliting the target data  using the covar table")

    #Get the typed chrom X
    X = config["input"]["x"]["typed"]

    #PLINK.py - Split the chromosome X per sex
    bfileFemale, bfileMale, X = separateGenotypedDataBySex(X, covarDict, args.folder, [],"BySex_ChrX", config["programs"]["plink1"], logFile, debug)

    #PCA.py - Calculate PCA and remove the outlier
    outlierList = detectPCAOutliersMaleFemale(bfileMale, bfileFemale, config, args.folder, args.threads, logFile, debug)

    #PLINK.py - Remove all outliers detected and build the Male and Female data
    nonOutlierFolder = f"{args.folder}/NonOutlierData"
    createFolder(nonOutlierFolder, log=False)
    bfileFemaleNonOut, bfileMaleNonOut, both = separateGenotypedDataBySex(X, covarDict, nonOutlierFolder, outlierList, "BySex_ChrX_nonOut",
                                                           config["programs"]["plink1"], logFile, debug)

    #PCA.py - Calculate PCA and remove the outlier for both (dataset with all males and female)
    outlierList = detectPCAOutliersBoth(both, config, args.folder, args.threads, logFile, debug)

    fileRemove = open(f"{args.folder}/OutlierDetection/removeBoth.txt", "w")
    for ID in outlierList:
        fileRemove.write(f"{ID}\t{ID}\n")
    fileRemove.close()

    #PLINK.py - Remove the outliers detected on the both dataset
    bfileBothNonOut = removeSamples(both, f"{args.folder}/OutlierDetection/removeBoth.txt", f"Both_ChrX_nonOut", nonOutlierFolder, config["programs"]["plink1"], logFile, debug, True)

    femalePC, malePC, bothPC = getRegressionPCA(bfileFemaleNonOut, bfileMaleNonOut, bfileBothNonOut, config, "AutoNonOutlier",args.folder, args.threads, logFile, debug)

    covarFileFemale, pheno, covarListFemale, updateFemale = buildCovarFileAndModel(femalePC, covarDict, config, nonOutlierFolder, "covarFemale", logFile, debug)
    covarFileMale, pheno, covarListMale, updateMale = buildCovarFileAndModel(malePC, covarDict, config, nonOutlierFolder, "covarMale", logFile, debug)
    covarFileBoth, pheno, covarListBoth, updateBoth = buildCovarFileAndModel(bothPC, covarDict, config, nonOutlierFolder, "covarBoth", logFile, debug)

    #Finally: Regression
    #First: prepare the pfile
    pfileMale = convertVCFToPFile(covarFileMale, config, nonOutlierFolder, "MaleData", updateMale, logFile, debug)
    pfileFemale = convertVCFToPFile(covarFileFemale, config, nonOutlierFolder, "FemaleData", updateFemale, logFile, debug)
    pfileBoth = convertVCFToPFile(covarFileBoth, config, nonOutlierFolder, "BothData",updateBoth, logFile, debug)

    rsquare = "0.5"
    firth = True
    regressionMale = runRegressionPlink2(pfileMale, covarFileMale, rsquare, firth, pheno, covarListMale, args.folder, "RegressionMale", config, args.threads, logFile, debug)
    regressionFemale = runRegressionPlink2(pfileFemale, covarFileFemale, rsquare, firth, pheno, covarListFemale, args.folder, "RegressionFemale", config, args.threads, logFile, debug)
    regressionBoth = runRegressionPlink2(pfileBoth, covarFileBoth, rsquare, firth, pheno, covarListBoth, args.folder, "RegressionBoth", config, args.threads, logFile, debug)

    gwamaMetaAnalysis(regressionFemale, regressionMale, pfileFemale, pfileMale, firth, "FemaleMale", f"{args.folder}/Results", config, logFile, debug)