import os, math, csv,sys
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
os.chdir("F:/bioinformatics/research_project/virus_sequences")

def blastAll(): #will search through subdirectories of a certain root and will send the fasta files to get blasted
    for path, subdirs, files in os.walk("F:/bioinformatics/research_project/virus_sequences"):#I plan on putting a prompt where you specify where your sequences are
        for name in files:
            if name[0:len(name)-8].find(".") == -1 :#so we don't analyse hidden files
                fileNames = os.path.join(path,name).replace('\\', '/')#clean the file name
                print(1)
                if fileNames.find(".fasta") != -1:#push fasta files only
                        print(blastn(fileNames))#this will blast your file
                        sys.exit()
    return["done blasting all"]

def blastn(name): #this function will takes a filepath of a fasta file as an argument and does a local blastn with the nt database
    outPutName = os.path.splitext(name)[0] #this is simply a variable to name the XML output file
    E_VALUE = pow(10,-40)
    blast = NcbiblastxCommandline(cmd = "blastn", query= name,
                                db="D:/database/nt/nt", evalue= E_VALUE,
                                outfmt=5, out= outPutName + ".xml",
                                maxtargetseqs= 100000)
    print('blasting '+ name)
    stdout,stderr = blast()

    return["done blasting " + name]


blastAll()
