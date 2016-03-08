#voici mon script pour trouver les termes GO de mes genes


import pymysql, os, uniprot, pprint, csv,urllib2, time,platform
from operator import itemgetter
from itertools import groupby

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")

#cette fonction va trouver mes fichier CSV qui contienne la liste des genes
#pour chacun des genomes de virus obtenu
def findCSV(): #will search through subdirectories of a certain root and will send the csv files
    for path, subdirs, files in os.walk("./research_project"):
        for name in files:
            if name[0:len(name)-4].find(".") == -1 :#so we don't analyse hidden files
                fileNames = os.path.join(path,name).replace('\\', '/')#clean the file name
                if fileNames.find("sequences.csv") != -1:#push the right csv file
                        outPutName = os.path.splitext(fileNames)[0]
                        print(outPutName+'.csv')#to follow along in the terminal
                        print(protIDret(outPutName))#s'occupe d'ouvrir le document
    return[]

#cette fonction fait des appelle mySQL au serveur de ebi afin d'obtenir les
#GO terms pour une gene. Celui-ci utilise les uniprotIDs pour le faire
def GetGoAnnotation(seqids):

    db = pymysql.connect(host = "mysql-amigo.ebi.ac.uk",
                     user = "go_select",
                     passwd = "amigo",
                     db = "go_latest",
                     port = 4085)

    cur = db.cursor()
    cur.execute(
          """
      SELECT
        term.name,
        term.acc,
        term.term_type
       FROM   gene_product
        INNER JOIN dbxref ON (gene_product.dbxref_id=dbxref.id)
        INNER JOIN species ON (gene_product.species_id=species.id)
        INNER JOIN association ON (gene_product.id=association.gene_product_id)
        INNER JOIN evidence ON (association.id=evidence.association_id)
        INNER JOIN term ON (association.term_id=term.id)
       WHERE
         dbxref.xref_key = %s;
         """, seqids)
    List = list()
    GO = list()
    f= cur.fetchall()
    for i in f:
        List.append(i[0] + ":" + i[2])
        GO.append(i[1])
    List = list(map(itemgetter(0), groupby(List)))
    GO = list(map(itemgetter(0), groupby(GO)))
    db.close()
    return[seqids,List,GO]

#cette fonction utilise un script ecrit par quelqu'un qui va prendre un
#Refseq_protein et va obtenir des uniprotID(si disponible)
def getUNIPid(line):
  seqid = list()
  seqid.append(line[11:line[0].find(']')-1])
  pairs = uniprot.batch_uniprot_id_mapping_pairs('P_REFSEQ_AC','ACC',seqid)
  reviewedpair = ""
  revStat= 0
  #Cette partie du code sert a determiner qu'elle des uniprotID sont reviewed
  #et va conserver le premier de la liste
  for y in pairs:
      uniCode = y[1]
      req = urllib2.Request('http://www.uniprot.org/uniprot/?query={}&sort=score&columns=reviewed&format=tab'.format(uniCode))
      web = urllib2.urlopen(req)
      for i in web:
          if i.find("Status") == -1:
              if i.find("unreviewed") == -1:
                  if revStat == 0 :
                      reviewedpair = y[1]
                      revStat = 1
                      break
              elif reviewedpair == "" :
                  reviewedpair = y[1]


  return[GetGoAnnotation(reviewedpair)]

#cette fonction ouvre les fichiers csv avec lesquels nous allons travaille
#il ouvre une version qui a ete GOed(to GO au passer) et le fichier d'interet
def protIDret(path):
    #ouverture du document que nous ecrivons
    with open(path +'_GOed.csv','wb') as p:
        writer = csv.writer(p)

        #ouverture du document qui contient les genes que nous voulons annote
        with open(path+'.csv', 'rb') as f:
            reader = csv.reader(f)
            your_list = list(reader)

        #ouverture de chaque element de la liste
        for i in  your_list:
            #ecriture du header
            if i[1].find("Accession") != -1:
                writer.writerow([ "Virus","Accession","Protein","Protein_id","Hits","UniprotID", "Terms","GO code"])
            else:
                #appelle a la fonction qui obtien les uniprotID pour chaque Refseq_protein
                for c,k,v in getUNIPid(i[3]):
                    i.append(c)
                    i.append(k)
                    i.append(v)
                writer.writerow(i)

    return[]



findCSV()
