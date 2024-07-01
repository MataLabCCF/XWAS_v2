import argparse
from Handlers.COVAR import readCovarFile
from Handlers.PLINK1 import *
from Utils.GENERAL import createFolder

# Differences from XWAS v1 and XWAS v2
# - We removed the country file. It was implemented on v1 and never used because the analysis per country
# requires a QC per country.


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PCA and regression')

    data = parser.add_argument_group("Data arguments")
    data.add_argument('-g', '--genotyped', help='Genotyped file name', required=False)
    data.add_argument('-i', '--imputed', help='Imputed file name', required=False)
    data.add_argument('-t', '--tableCovar', help='File with covariatives to be added to the model', required=True)
    data.add_argument('-r', '--reference', help='Data with parental reference data to run PCA (optional)',
                      required=False, default="")

    output = parser.add_argument_group("Output arguments")
    output.add_argument('-n', '--name', help='Analysis name', required=True)
    output.add_argument('-f', '--folder', help='Folder to output files', required=True)

    programs = parser.add_argument_group("Programs arguments")
    programs.add_argument('-G', '--gwama', help='GWAMA program (default = gwama)', required=False, default="gwama")
    programs.add_argument('-p', '--plink2', help='Path of PLINK 2 (default = plink2)', required=False, default="plink2")
    programs.add_argument('-k', '--plink1', help='Path of PLINK 1 (default = plink)', required=False, default="plink")

    programs.add_argument('-P', '--python', help='Path of Python 3 (default = python)', required=False,
                          default="python")
    programs.add_argument('-R', '--runPCA', help='Path of runPCA script (default = runPCA.R)', required=False,
                          default="runPCA.R")

    args = parser.parse_args()

    logFile = createFolder(args.folder)

    covarDict = readCovarFile(args.tableCovar)
    bfileFemale, bfileMale = separateGenotypedDataBySex(args.genotyped, covarDict, args.name, args.folder, args.plink1, logFile)

    if args.reference != "":
        bfileMaleRef, bfileFemaleRef = separateGenotypedDataBySexWithoutCovar(args.genotyped, "Ref", args.folder, args.plink1, logFile)

        bfileFemale = mergeRefAndTarget(bfileFemale, bfileFemaleRef, f"{args.name}_mergedWithRef_Female", args.plink1, logFile)
        bfileMale = mergeRefAndTarget(bfileMale, bfileMaleRef, f"{args.name}_mergedWithRef_Male", args.plink1, logFile)