
import os,platform, subprocess, sys, Levenshtein, csv, multiprocessing, os.path
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")

global infoList
infoList = []

def compareAlign(fileName):

    splitName = fileName.split("/")

    geneList = [splitName[4],splitName[6][:splitName[6].find("_")],
                splitName[7]]

    foldName = fileName[:fileName.rfind('/')+1]
    for file in os.listdir(foldName):
        if '_refseq.fasta' in file:
            ref = foldName+file

    out = fileName[:fileName.rfind("_")] + "_filtered.fasta"
    alignEm(ref,fileName,out)
    if os.path.isfile(out):
        with open(out, "r") as fasta:
            sequences = fasta.read()
            freq = sequences.count(">")
    else: freq = 0
    
    geneList.append(freq)
    infoList.append(geneList)

    return

def alignEm(refSeq, record, out):
    shell = "F:/bioinformatics/usearch -usearch_global " + record + " -db " + refSeq + " -strand plus -id 0.97 -matched " + out
    subprocess.call(shell,shell=True)

    


def crawl(folder):
    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')   
            if "_Fetch.fasta" in fileName:

                compareAlign(fileName)

    return

crawl("./research_project/test")
with open("./research_project/test/filtered_seq_info_test.csv","wb") as f:
    writer = csv.writer(f)
    writer.writerow(["type","virus","gene","number sequence"])
    for i in infoList:
        writer.writerow(i)

