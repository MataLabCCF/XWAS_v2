import os
import sys
sys.path.append('../')
from Handlers.PLINK1 import countFam, countBim
from Handlers.PLINK2 import countPsam, countPvar

def createFolder(folder, log=True):
    if not os.path.exists(folder):
        os.mkdir(folder)
    if log:
        logFile = open(f"{folder}/XWAS.log", "w")
        return logFile

def execute(commandLine, type, outputPrefix, logFile):
    os.system(commandLine)

    writeInfoFile = True

    if type == "bfile":
        numSample = countFam(outputPrefix)
        numVar = countBim(outputPrefix)
    elif type == "pfile":
        numSample = countPsam(outputPrefix)
        numVar = countPvar(outputPrefix)
    else:
        writeInfoFile = False

    logFile.write("Command line:\n")
    logFile.write(f"\t{commandLine}:\n")
    logFile.write(f"\n")
    if writeInfoFile:
        logFile.write(f"Output information:\n")
        logFile.write(f"\tNumber of samples: {numSample}\n")
        logFile.write(f"\tNumber of varaintss: {numVar}\n")
        logFile.write(f"\n")
    logFile.write("=================================================================================================\n")

