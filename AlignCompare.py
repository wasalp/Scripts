
import os,platform, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")


def compareAlign(fileName):
    FastaRec = []
    foldName = fileName[:fileName.rfind('/')+1]
    for file in os.listdir(foldName):
        if '_refseq.fasta' in file:
            ref = foldName+file
    refSeq = SeqIO.parse(ref,"fasta").next()
    record = SeqIO.parse(fileName,"fasta")
    print fileName
    alignEm(refSeq,record)
    return

def alignEm(refSeq, record):

    for i in record:
        for a in pairwise2.align.globalxx(refSeq.seq, i.seq):
            print (format_alignment(*a))
        break


    return


def crawl(folder):

    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if ".fasta" in fileName and ".DS_Store" not in fileName:
                if "_genes" in fileName and ".gb" not in fileName:
                    if "_refseq" not in fileName:
                        compareAlign(fileName)
    return


crawl("./pretest_folder")
