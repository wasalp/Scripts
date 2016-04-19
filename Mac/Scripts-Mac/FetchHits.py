import sys,os,csv,time, platform
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import Entrez

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")
def GetFiles():

    for path, subdirs, files in os.walk("./project_temp"):#I plan on putting a prompt where you specify where your sequences are
        for name in files:
            if name[0:len(name)-8].find(".") == -1 :#so we don't analyse hidden files
                PathName = path.replace('\\','/')
                fileNames = os.path.join(path,name).replace('\\', '/')#clean the file name
                if fileNames.find(".xml") != -1:#I think I can change this to if "XML" in fileNames:
                        outPutName = os.path.splitext(fileNames)[0] #this is simply a variable to name the XML output file
                        print(outPutName+'.xml')
                        XMLparse(PathName, outPutName+'.xml')
    return

def XMLparse(Path,File):
    E_VALUE = pow(10,-20)
    result = NCBIXML.parse(open(File))
    #gi|353260843|gb|JN555585.1|
    r = os.path.splitext(File)[0]
    s = r[r.rfind('/'):]
    print s
    Path = Path + s[: s.find('_')+1] + 'genes'
    try:
        os.mkdir(Path)
    except:
        print "already exists"

    for blast_record in result:
#create a folder for all the genes of our virus(probably at a previous point in the script)
        queryName = blast_record.query
        protID = queryName[queryName.find("protein_id="):queryName.find("]", queryName.find("protein_id="))]
        FoldName =  protID[11:]
        print FoldName
        print Path
        if ("NP_062890.1" in FoldName or "NP_062888.1" in FoldName or "NP_056796.1" in FoldName
          or "NP_056794.1" in FoldName or "NP_066246.1" in FoldName or "NP_066244.1" in FoldName
          or "NP_041332.1" in FoldName or "NP_041327.1" in FoldName or "YP_003708381.1" in FoldName
          or "YP_003708382.1" in FoldName or "YP_001129462.1" in FoldName or "YP_001129465.1" in FoldName):
            try:

    #create a folder for each individual gene(extract either from each blast record(if possible) or from the fasta file
    #get multi-genbank file of each hit with the accession number(already works)
            #List of accession number from BLAST hits
                AccList = list()
                for alignment in blast_record.alignments:
                    AccList.append(alignment.hit_id.split("|")[3])
                os.mkdir(Path + '/' + FoldName)
                handlez = FetchGB(AccList)

                with open(Path + '/' + FoldName + '/' + FoldName + "_Fetch.gb",'wb') as f:

                    for line in handlez:
                        f.write(line)
                time.sleep(15)
            except:
                if len(FoldName) < 100:
                    if os.path.isfile(Path + '/' + FoldName + '/' + FoldName + "_Fetch.gb"):
                        if os.stat(Path + '/' + FoldName + '/' + FoldName + "_Fetch.gb").st_size < 300:

                            handlez = FetchGB(AccList)
                            with open(Path + '/' + FoldName + '/' + FoldName + "_Fetch.gb",'wb') as f:

                                for line in handlez:
                                    f.write(line)
                    else:
                        handlez = FetchGB(AccList)
                        with open(Path + '/' + FoldName + '/' + FoldName + "_Fetch.gb",'wb') as f:
                            for line in handlez:
                                f.write(line)
                print "already exists"
        #CDSEx(Path+ '/' + FoldName + "_Feth.gb", protein[7:])
    result.close()
    return

def CDSEx(NameFile, protein):
    #for seq_record in SeqIO.parse(NameFile, "genbank"):
    print "1"


    return []

def FetchGB(AccList):

    db = "nucleotide"
    Entrez.email = "louisparent@gmail.com"
    batchSize = 3000
    retmax = 10**9



    query = " ".join(AccList)
    handle = Entrez.esearch( db=db,term=query,retmax=retmax )
    while True:
        try:
            print query
            if not query:
                break
            sys.stderr.write( "Fetching entries from GenBank")

            giList = Entrez.read(handle)['IdList']
            sys.stderr.write( "Found %s GI: %s\n" % (len(giList), ", ".join(giList[:10])))
    #post NCBI query
            search_handle     = Entrez.epost(db=db, id=",".join(giList))
            search_results    = Entrez.read(search_handle)
            webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"]
    #fecth all results in batch of batchSize entries at once
            for start in range( 0,len(giList),batchSize ):
                  sys.stderr.write( " %9i" % (start+1,))
      #fetch entries in batch
                  handle = Entrez.efetch(db=db, rettype="gb", retstart=start,
                                      retmax=batchSize, webenv=webenv, query_key=query_key)

      #print output to stdout
            break
        except:
            print 'try again'
            time.sleep(2)
    return handle.read()

GetFiles()
