import os, subprocess, time
import Tkinter as tk
from tkFileDialog import askdirectory

os.chdir("/Users/louis/Desktop/exempleMEME")
def treethis(fileName,analysis):
    print("queuing")
    outName = getout(fileName,analysis)
    treepath = crawl(fileName[:fileName.rfind("/")+1],".tree",init=1)
    cmd = "HYPHYMP CPU=4 BASEPATH=/Users/louis/Desktop/exempleMEME/ QuickSelectionDetection_example.bf"
    cmd += "<<<$'" + fileName + "\n" + treepath + "\n" +outName+"\n"+analysis +"'&"
    subprocess.call(cmd,shell=True)

def timekeeper(fileName,analysis):
    while True:
        if int(subprocess.check_output("ps aux|grep HYPHY|wc -l",shell=True)) > 4:
            time.sleep(120)
        else:
            treethis(fileName,analysis)
            break

    return



def getout(fileName,analysis):
    if not os.path.exists(fileName[:fileName.rfind("/")+1] + analysis):
        os.makedirs(fileName[:fileName.rfind("/")+1] + analysis)
    splitName = fileName.split("/")
    outName = (fileName[:fileName.rfind("/")+1] + analysis + "/" + splitName[7])
    return outName

def crawl(folder, whitelist,blacklist="aaaaaaaaaaaaaaaaaaa", init=0,analysis=""):

    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if whitelist in fileName and blacklist not in fileName:
                if init == 0:
                    timekeeper(fileName,analysis)
                elif init == 1:
                    return fileName

root = tk.Tk()
root.withdraw()
directory = askdirectory()
root.destroy()

crawl(directory, "rn.fasta.nt_cleanali_rn.fasta", ".tree", analysis="MEME")
crawl(directory, "rn.fasta.nt_cleanali_rn.fasta", ".tree", analysis="SPIDER")
