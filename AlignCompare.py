
import os,platform, subprocess, sys, Levenshtein, csv, multiprocessing
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")


def compareAlign(fileName):
    print fileName
    with open(os.path.splitext(fileName)[0] + "_filtered.fasta",'wb') as h:
        FastaRec = []
        foldName = fileName[:fileName.rfind('/')+1]
        for file in os.listdir(foldName):
            if '_refseq.fasta' in file:
                ref = foldName+file
        refSeq = SeqIO.parse(ref,"fasta").next()
        record = SeqIO.parse(fileName,"fasta")
        score =[]

        fileName = foldName + fileName[fileName.rfind('/')+1:]
        for i in record:
            if alignEm(refSeq,i) >= 0.7:
                print i
                FastaRec.append(i)

        SeqIO.write(FastaRec, os.path.splitext(fileName)[0] + "_filtered.fasta" ,"fasta")
    return

def alignEm(refSeq, record):
    try:
        muscle_cline = MuscleCommandline("../muscle",maxiters= 1, diags=True)
        child = subprocess.Popen(str(muscle_cline),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=False,
                                shell=True)

        SeqIO.write([refSeq,record],child.stdin,"fasta")
        child.stdin.close()
        align = AlignIO.read(child.stdout,"fasta")
        return Levenshtein.ratio(str(align[0].seq), str(align[1].seq))
    except:
        return 0


def crawl(folder):
    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if ".fasta" in fileName and ".DS_Store" not in fileName:
                if "_genes" in fileName and ".gb" not in fileName and "_filtered" not in fileName:
                    if fileName.find("refseq") == -1 :
                        compareAlign(fileName)

    return


crawl("./research_project")
