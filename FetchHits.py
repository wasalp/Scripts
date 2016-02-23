import sys,os,csv,time
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import Entrez
os.chdir("F:/bioinformatics")

def GetFiles():

    for path, subdirs, files in os.walk("F:/bioinformatics/research project/virus sequences"):#I plan on putting a prompt where you specify where your sequences are
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
    E_VALUE = pow(10,-5)
    result = NCBIXML.parse(open(File))
    #gi|353260843|gb|JN555585.1|
    r = os.path.splitext(File)[0]
    s = r[r.rfind('/'):]
    print s
    Path = Path + s[: s.find('_')+1] + 'genes'
    try:
        os.mkdir(Path)
    except WindowsError:
        print "already exists"
        
    for blast_record in result:
#create a folder for all the genes of our virus(probably at a previous point in the script)
        queryName = blast_record.query
        protein = queryName[queryName.find("protein=")-1:queryName.find("]", queryName.find("protein="))+1]
        protein = protein.replace('/', '_')
        protID = queryName[queryName.find("protein_id=")-1:queryName.find("]", queryName.find("protein_id="))+1]
        if len(protein) > 150:
            protein = queryName[queryName.find("gene=")-1:queryName.find("]", queryName.find("gene="))+1]
        FoldName = '[' + protein[protein.find("=")+1:] + '[' + protID[12:]
        print FoldName
        print Path
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
        except WindowsError:
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
    batchSize = 20
    retmax = 10**9



    query = " ".join(AccList)
    while True:
        try:
            sys.stderr.write( "Fetching entries from GenBank")
            handle = Entrez.esearch( db=db,term=query,retmax=retmax )
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
