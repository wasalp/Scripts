#The purpose of this function will be to compare a series of coding sequences
#that was fetched from NCBI databases to a reference sequence with the goal of
#of putting all similar sequences(with a certain identity treshold) in the same
#multi-fasta file for further analysis.

import os,platform
from Bio import SeqIO
import Bio.Align.Applications
import subprocess

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")

#os.chdir("./pretest folder")

#this function will get the reference sequence
def getRefSeq(fileName, identifier):
    geneFolder = fileName[:fileName.rfind('/')+1] + identifier + "_genes"
    if os.path.isdir(geneFolder):
        record = SeqIO.parse(fileName,'fasta')
        dirList = os.listdir(geneFolder)
        for seq_record in record:
            ID = fetchID(seq_record.description)
            for directory in dirList:
                if ID in directory:
                    fastaFile = (geneFolder + '/' + directory + '/' + ID
                                + "_refseq.fasta")
                    SeqIO.write(seq_record, fastaFile, "fasta")

    return



def crawl(folder):

    for path,subdirs,files in os.walk(folder):
        for name in files:
            if name[0:len(name)-8].find(".") == -1:
                fileName = os.path.join(path,name).replace('\\', '/')
                if ".fasta" in fileName: getRefSeq(fileName,fetchIden(fileName))
    return

def fetchID(ID):

    ID = ID[ID.find("protein_id")+11:]
    ID = ID[:ID.find(']')]

    return ID

def fetchIden(fileName):

    identifier = os.path.splitext(fileName)[0]
    identifier = identifier[identifier.rfind('/')+1:]
    identifier = identifier[:identifier.find('_')]

    return identifier

#Hopefully this function will be able to extract cds features from genbank files
#using the readseq java program
def readseqGB():





  return


print os.getcwd()
crawl("./pretest_folder")
