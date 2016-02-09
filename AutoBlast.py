import os, math, csv
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
os.chdir("D:/bioinformatics/research project/virus sequences")

def blastAll(): #will search through subdirectories of a certain root and will send the fasta files to get blasted
    for path, subdirs, files in os.walk("D:/bioinformatics/research project"):#I plan on putting a prompt where you specify where your sequences are
        for name in files:
            if name[0:len(name)-8].find(".") == -1 :#so we don't analyse hidden files
                fileNames = os.path.join(path,name).replace('\\', '/')#clean the file name
                if fileNames.find(".fasta") != -1:#push fasta files only
                        print(blastn(fileNames))#this will blast your file
    return["done blasting all"]

def blastn(name): #this function will takes a filepath of a fasta file as an argument and does a local blastn with the nt database
    outPutName = os.path.splitext(name)[0] #this is simply a variable to name the XML output file
    blast = NcbiblastxCommandline(cmd = "blastn", query= name, db="D:/database/nt/nt", evalue=0.2, outfmt=5, out= outPutName + ".xml")
    print('blasting')
    stdout,stderr = blast()

    return["done blasting " + name]


blastAll()
