#simple script to extract all CDS from information from genbank fileNames

import os,platform, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")

readseqLoc = "/Users/louis/Desktop/readseq.jar"

def getGbInf(fileName):
    try:
        record = SeqIO.parse(fileName,"genbank")
        FastaRec = []
        for i in record:
            if i.features:
                for feature in i.features:
                    if feature.type == "CDS":
                        FastaRec.append(
                            SeqRecord(feature.location.extract(i).seq,
                            id=str(feature.qualifiers.get("protein_id", "???"))))

        SeqIO.write(FastaRec, os.path.splitext(fileName)[0]+".fasta","fasta")
    except:
        print fileName


    return



def crawl(folder):

    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if ".gb" in fileName and ".DS_Store" not in fileName:
                getGbInf(fileName)
    return

crawl("./research_project")
