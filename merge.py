import csv,os

os.chdir("F:/bioinformatics")
global newlist
newlist = []


def second(line):
    line[3] = line[3][line[3].find("=")+1:]
    line[4] = line[4][line[4].find("=")+1:]
    with open("./research_project/test/filtered_seq_info_test.csv", "rb") as p:
        seqinfo = csv.reader(p)
        for row in seqinfo:
            if row[2] == line[4][line[4].find("=")+1:]:
                line.append(row[3])
                newlist.append(line)
    
   

with open("./research_project/allMerged.csv", "rb") as f:
        allmerged = csv.reader(f)
        allmerged.next
        newlist= [["Type","Virus","Accession","Protein","Protein_id",
                        "Hits","uniprotID","Terms","GO code","seqnum"]]
        for lines in allmerged:
            second(lines)

with open("./research_project/info.csv","wb") as f:
    writer = csv.writer(f)
    for i in newlist:
        writer.writerow(i)
