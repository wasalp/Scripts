from Tkinter import Tk
from tkFileDialog import askdirectory
import csv,os
import os,platform
from Bio import SeqIO



def extract(fileName,tresholdSPIDER,tresholdMEME):

    Dict = {}
    Matchlist = {}
    splitFileName = fileName.split("/")
    with open(fileName,"rb") as f:
        reader = csv.reader(f)
        reader.next()
        for lines in reader:
            if lines[0] not in Dict:
                if float(lines[4]) > tresholdSPIDER:
                    Dict[lines[0]] = float(lines[4])


    with open(searchDir(fileName[:fileName.rfind("/")+1],".branches"),"rb") as p:
        reader = csv.reader(p)
        reader.next()
        reader.next()
        reader.next()
        for lines in reader:
            if lines[0] not in Matchlist:
                if float(lines[2]) > tresholdMEME:
                    Matchlist[lines[0]] = lines[2]



    record = SeqIO.parse(searchDir(fileName[:fileName.rfind("/")+1],".fasta.nt_cleanali.fasta"),"fasta")
    seqlen = len(record.next())
    SeqNum = 1
    for i in record:
        SeqNum += 1

    matchcount = 0
    for i in Matchlist:
        if i in Dict:
            matchcount += 1

    Alllist.append([splitFileName[5],splitFileName[6],
                    splitFileName[7],splitFileName[8],len(Dict),len(Matchlist),
                    matchcount,tresholdSPIDER,tresholdMEME])



    return


def searchDir(folder, keyword):
    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if keyword in fileName:
                return fileName

def crawl(folder, keyword,tresholdSPIDER,tresholdMEME):
    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if keyword in fileName:
                extract(fileName,tresholdSPIDER,tresholdMEME)

Tk().withdraw()
global directory
directory = askdirectory()
global Alllist
Alllist = []
for y in range(10,20,1):
    b = float(y)
    for i in range(10,20,1):
        j = float(i)
        crawl(directory,"Edge_support.csv",j/20,b/20)


with open(directory+"/AllCoEvSel.csv","wb") as allCo:
    allCoWriter = csv.writer(allCo)
    allCoWriter.writerow(["type","virus","gene","function","SPIDER","MEME","intersection","tresholdSPIDER","tresholdMEME"])
    for i in Alllist:
        allCoWriter.writerow(i)
