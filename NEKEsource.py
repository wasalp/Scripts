#simple script to extract all CDS from information from genbank fileNames

import os,platform, subprocess,csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")


def getGbInf(fileName):
    with open("Neke.csv", "wb") as f:
        writer = csv.writer(f)
        record = SeqIO.parse(fileName,"genbank")
        writer.writerow(["Accession","collection_date"])
        for i in record:
            if i.features:
                for feature in i.features:
                    if feature.type == "source":
                        temp = str(feature.qualifiers['collection_date']).replace("['", "")
                        temp = temp.replace("']","")
                        writer.writerow([i.id,temp])


    return



def crawl(folder):

    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if ".gb" in fileName and ".DS_Store" not in fileName:
                getGbInf(fileName)
    return

crawl("./")
