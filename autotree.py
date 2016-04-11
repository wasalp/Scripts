import os, subprocess
from Tkinter import Tk
from tkFileDialog import askdirectory


def treethis(fileName):

    outName = getout(fileName)
    cmd = "~/Desktop/FastTree -gtr -nt -nosupport " + fileName + ">" + outName
    print(cmd)
    subprocess.call(cmd,shell=True)


def getout(fileName):

    splitName = fileName.split("/")
    outName = (fileName[:fileName.rfind("/")+1] + splitName[5] + ".tree")
    return outName

def crawl(folder, whitelist,blacklist):
    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if whitelist in fileName and blacklist not in fileName:
                treethis(fileName)


Tk().withdraw()
directory = askdirectory()

crawl(directory, "rn.fasta.nt_cleanali.fasta", ".tree")
