def countPsam(filePrefix):
    famFile = open(f"{filePrefix}.psam")
    count = 0

    header = True
    for line in famFile:
        if header:
            header = False
        else:
            count = count + 1
    famFile.close()
    return count

def countPvar(filePrefix):
    bimFile = open(f"{filePrefix}.bim")
    count = 0

    header = True
    for line in bimFile:
        if header:
            if "#CHROM" in line:
                header = False
        else:
            count = count + 1

    bimFile.close()
    return count
