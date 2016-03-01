
import os,platform, subprocess, sys, Levenshtein, csv
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline

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
    score =[]
    x=0
    fileName = fileName.replace("][",'_')
    fileName = fileName.replace('[','')
    fileName = fileName.replace(']','')
    with open(os.path.splitext(fileName)[0]+'.csv','wb') as f:
        writer = csv.writer(f)
        writer.writerow(["index","score"])
        for i in record:
            simi = [x , alignEm(refSeq,i)]
            writer.writerow(simi)

    return

def alignEm(refSeq, record):
    try:
        muscle_cline = MuscleCommandline("../muscle",maxiters= 1, diags=True)
        child = subprocess.Popen(str(muscle_cline),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True,
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
                if "_genes" in fileName and ".gb" not in fileName:
                    if "_refseq" not in fileName:
                        compareAlign(fileName)
    return


crawl("./research_project")
