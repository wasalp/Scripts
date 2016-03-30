from Tkinter import Tk
from tkFileDialog import askdirectory
import csv,os





def extract(fileName):
    Dict = {}
    Plist = []
    Dlist = []
    with open(fileName,"rb") as f:
        reader = csv.reader(f)
        reader.next()
        reader.next()
        reader.next()
        for lines in reader:
            if lines[0] in Dict:
                Plist.append(lines[2] +":"+lines[3])
                Dict[lines[0]] = Plist
            else:
                Plist = []
                Plist.append(lines[2] +":"+lines[3])
                Dict[lines[0]] = Plist
    with open(fileName[:fileName.rfind("/")+1]+"Edge_support.csv") as f:
        reader = csv.reader(f)
        reader.next()
        for lines in reader:
            if lines[0] in Dict:
                Dlist.append(lines[2] +":"+lines[3])
                Dict[lines[0]] = Plist
            else:
                Plist = []
                Plist.append(lines[4])
                Dict[lines[0]] = Plist
                
        for lines in reader
        for i in range(5):
            dictlist = Dict.items()
            print(dictlist[i])
            print("/n")



def crawl(folder):
    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            if ".branches" in fileName:
                extract(fileName)

    return




Tk().withdraw()
directory = askdirectory()
crawl(directory)
