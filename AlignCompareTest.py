
import os,platform, subprocess, sys, Levenshtein, csv
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastxCommandline

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")


def compareAlign(fileName):
    FastaRec = []
    x=0
    foldName = fileName[:fileName.rfind('/')+1]
    for file in os.listdir(foldName):
        if '_refseq.fasta' in file:
            ref = foldName+file
    refSeq = SeqIO.parse(ref,"fasta").next()
    record = SeqIO.parse(fileName,"fasta")
    score =[]

    for i in record:
        x = x+1
        alignEm(refSeq,i)
    print x
    sys.exit()
    return

def alignEm(refSeq, record):
    blast = NcbiblastxCommandline(cmd = "blastn", outfmt=5)
    child = subprocess.Popen(str(blast),
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    SeqIO.write(record,child.stdin,"fasta")
    out,err = process.communicate(subject = refSeq)
    print out
    print err
    sys.exit()

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
