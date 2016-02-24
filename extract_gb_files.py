#simple script to extract all CDS from information from genbank fileNames

import os,platform, subprocess

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")

readseqLoc = "/Users/louis/Desktop/readseq.jar"

def getGbInf(fileName):
    os.chdir(fileName[:fileName.rfind('/')+1]
    callText = "java -cp {0} run -f=8  -o= {1} -feat=CDS {2}".format(
        readseqLoc,os.path.splitext(fileName)[0]+".fasta",fileName)
    print callText
    subprocess.call(callText,shell=True)

    return



def crawl(folder):

    for path,subdirs,files in os.walk(folder):
        for name in files:
            fileName = os.path.join(path,name).replace('\\', '/')
            fileName = fileName.replace("]","\]")
            fileName = fileName.replace('[','\[')
            fileName = fileName.replace(' ', '\ ')
            if ".gb" in fileName and ".DS_Store" not in fileName:
                getGbInf(fileName)
    return

crawl("./pretest folder")
