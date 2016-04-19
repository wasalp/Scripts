
import os,platform, subprocess, sys, Levenshtein, csv, multiprocessing
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
import Tkinter as tk
from tkFileDialog import askdirectory

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")


def reName(fileName):
    print fileName
    with open(fileName,"r") as fasta:
        with open(os.path.splitext(fileName)[0]+"_rn.fasta","wb") as new:
            for line in fasta:
                if line.startswith('>'):
                    line = line.replace("['","")
                    if "]" in line:
                        line = line[:line.find("]")-1]
                    if "." in line:
                        line = line[:line.find(".")+1]
                    line = line[:line.rfind("/")]
                    line = line.replace("'","")
                    print line
                    new.write(line+"\n")
                else:
                    new.write(line)

    return


def crawl(folder, whitelist,blacklist="aaaaaaaaaaaaaaaaaaa"):
    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if whitelist in fileName and blacklist not in fileName:
                reName(fileName)

    return

root = tk.Tk()
root.withdraw()
directory = askdirectory()
root.destroy()

crawl(directory,"fasta.nt_cleanali.fasta")
