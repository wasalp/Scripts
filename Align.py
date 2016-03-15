
import os,platform, subprocess, sys, Levenshtein, csv, multiprocessing
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")


def reName(fileName):
    print fileName
    with open(fileName,"r") as fasta:
        with open(os.path.splitext(fileName)[0]+"_clean.alignment","wb") as new:
            for line in fasta:
                if line.startswith('>'):
                    line = line.replace("['","")
                    line = line[:line.find("]")-1]
                    line = line[:line.find(".")]
                    print line
                    new.write(line+ "\n")
                else:
                    new.write(line)

    return


def crawl(folder):
    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if ".alignment" in fileName and ".DS_Store" not in fileName and "clean" not in fileName:
                reName(fileName)

    return


crawl("./project_temp")