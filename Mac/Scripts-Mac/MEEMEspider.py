from Tkinter import Tk
from tkFileDialog import askdirectory
import csv,os
import os,platform
from Bio import SeqIO



def extract(fileName,treshold):

    Dict = {}
    Matchlist = {}
    Plist = []
    splitFileName = fileName.split("/")
    with open(fileName,"rb") as f:
        reader = csv.reader(f)
        reader.next()
        for lines in reader:
            if lines[0] in Dict:
                Plist.append(float(lines[4]))
                Dict[lines[0]] = Plist
            else:
                Plist = []
                Plist.append(float(lines[4]))
                Dict[lines[0]] = Plist


    with open(searchDir(fileName[:fileName.rfind("/")+1],".branches"),"rb") as p:
        reader = csv.reader(p)
        reader.next()
        reader.next()
        reader.next()
        for lines in reader:
            Matchlist[lines[0]] = 0

    record = SeqIO.parse(searchDir(fileName[:fileName.rfind("/")+1],".fasta.nt_cleanali.fasta"),"fasta")
    seqlen = len(record.next())
    SeqNum = 1
    for i in record:
        SeqNum += 1
    print(SeqNum)

    matchcount = 0
    with open(fileName[:fileName.rfind("/")+1]+ "CoEvSel.csv", "wb") as f:
        writer = csv.writer(f)
        for i in Matchlist:
            if i in Dict:
                for vals in Dict[i]:
                    if vals > treshold :
                        matchcount += 1

    with open(fileName[:fileName.rfind("/")+1]+ "CoEvSel.csv", "wb") as f:
        writer = csv.writer(f)
        for i in Matchlist:
            if i in Dict:
                for vals in Dict[i]:
                    if vals > treshold :
                        print([i,vals])
                        writer.writerow([i,vals])
                        Alllist.append([splitFileName[5],splitFileName[6],
                                 splitFileName[7],i,vals,seqlen,SeqNum,matchcount])


    return


def searchDir(folder, keyword):
    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if keyword in fileName:
                return fileName

def crawl(folder, keyword):
    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if keyword in fileName:
                extract(fileName,0.80)

Tk().withdraw()
global directory
directory = askdirectory()
global Alllist
Alllist = []
crawl(directory,"Edge_support.csv")
print(Alllist)
with open(directory+"/AllCoEvSel.csv","wb") as allCo:
    allCoWriter = csv.writer(allCo)
    allCoWriter.writerow(["type","virus","gene","site","value","seqlen","seqnum","hitnum"])
    for i in Alllist:
        allCoWriter.writerow(i)
