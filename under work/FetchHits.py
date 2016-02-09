import sys
from Bio import Entrez

db = "nuccore"
Entrez.email = "louisparent@gmail.com"
batchSize = 100
retmax = 10**9

query = " ".join(["KP965141", "KP965140"])
print query

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
  handle = Entrez.efetch(db=db, rettype="gb", retstart=start, retmax=batchSize, webenv=webenv, query_key=query_key)
  #print output to stdout
  sys.stdout.write(handle.read())
