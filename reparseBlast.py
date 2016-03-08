import os, math, csv,platform
from Bio.Blast import NCBIXML

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")

def startReparse(): #will search through subdirectories of a certain root and will send the fasta files to get blasted
    for path, subdirs, files in os.walk("./research_project"):#I plan on putting a prompt where you specify where your sequences are
        for name in files:
            if name[0:len(name)-8].find(".") == -1 :#so we don't analyse hidden files
                fileNames = os.path.join(path,name).replace('\\', '/')#clean the file name
                if fileNames.find(".xml") != -1:#push fasta files only
                        outPutName = os.path.splitext(fileNames)[0] #this is simply a variable to name the XML output file
                        print(outPutName+'.xml')
                        XMLparse(outPutName+'.xml')
    return[]


def XMLparse(Path):
    E_VALUE = pow(10,-20)
    result = NCBIXML.parse(open(Path))
    numBlastHits = 0
    hitsDict = {}
    for blast_record in result:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect <= E_VALUE:
                    numBlastHits += 1
        hitsDict[blast_record.query] = numBlastHits
        numBlastHits = 0
    result.close()
    print(writeCSV(Path, hitsDict))
    return[]

def fetchIden(fileName):

    identifier = os.path.splitext(fileName)[0]
    identifier = identifier[identifier.rfind('/')+1:]
    identifier = identifier[:identifier.find('_')]

    return identifier

def writeCSV(Path, Dict):
    CleanKey = ""
    CleanerKey = ''
    RowList = []
    writer = csv.writer(open(os.path.splitext(Path)[0]+'.csv','wb'))
    writer.writerow(["Virus","Accession","Protein","Protein_id", "BlastHits", "uniprotID"])
    for key, value in Dict.items():
        if key.find('lcl') != -1 :
            genePos = key.find('protein=')
            protein_id = key.find('protein_id=')
            RowList = [ fetchIden(Path),
                        str( key[4:key.find('_',key.find('_')+1)] ),
                        str( key[genePos:key.find(']',genePos)] ),
                        str( key[protein_id:key.find(']',protein_id)] )
                        ,value ]
        else:
            genePos = key.find("product=")
            protein_id = key.find("protein_id")
            RowList = [ fetchIden(Path),
                        str( key[0:key.find('_',key.find('_')+1)] ),
                        str(  key[genePos:key.find(';',genePos)] ),
                        str( key[protein_id:key.find(';',protein_id)] )
                        ,value]
            for i in range(len(RowList)-1):
                CleanerKey = RowList[i].replace('"','')
                CleanKey = CleanerKey.replace("'","")
                RowList[i] = CleanKey
        writer.writerow(RowList)

    return["done"]

startReparse()
