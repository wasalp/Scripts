#This is a simple program that will find all csv files in subfolders and append
#them together

import csv,os


os.chdir("/Users/louis/Desktop/bioinformatics/research project")


def protIDret():
    with open('/Users/louis/Desktop/bioinformatics/research project/allMerged.csv','wb') as p:
        writer = csv.writer(p)
        writer.writerow(["Type","Accession","Protein","Protein_id","Hits","uniprotID","Terms","GO code"])
        for path, subdirs, files in os.walk("/Users/louis/Desktop/bioinformatics/research project/virus sequences"):#I plan on putting a prompt where you specify where your sequences are
            for name in files:
                if name[0:len(name)-4].find(".") == -1 :#so we don't analyse hidden files
                    fileNames = os.path.join(path,name).replace('\\', '/')#clean the file name
                    if fileNames.find("GOed.csv") != -1:
                        #if fileNames.find(".csv") != -1:#push csv
                            outPutName = os.path.splitext(fileNames)[0] #this is simply a variable to name the csv output file
                            if outPutName.find("Users") != -1:
                                ViralType = outPutName.split("/")[7]
                            else:
                                ViralType = outPutName.split("/")[4]
                            print(ViralType)
                            print(outPutName+'.csv')#to follow along in the terminal


                            with open(outPutName +'.csv', 'rb') as f:
                                reader = csv.reader(f)
                                your_list = list(reader)

                            for i in  your_list:
                                if i[0].find("Accession") == -1:
                                    i.insert( 0 , ViralType )
                                    writer.writerow(i)

    return[]



protIDret()
