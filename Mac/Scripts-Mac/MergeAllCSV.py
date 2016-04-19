#This is a simple program that will find all csv files in subfolders and append
#them together

import csv,os,platform

if 'Darwin' in platform.system():
    os.chdir("/Users/louis/Desktop/bioinformatics/")
elif 'Windows' in platform.system():
    os.chdir("F:/bioinformatics")


def protIDret():
    with open('./research_project/allMerged.csv','wb') as p:
        writer = csv.writer(p)
        writer.writerow(["Type","Virus","Accession","Protein","Protein_id","Hits","uniprotID","Terms","GO code"])
        for path, subdirs, files in os.walk("./research_project/virus_sequences"):#I plan on putting a prompt where you specify where your sequences are
            for name in files:
                if name[0:len(name)-4].find(".") == -1 :#so we don't analyse hidden files
                    fileNames = os.path.join(path,name).replace('\\', '/')#clean the file name
                    if fileNames.find("GOed.csv") != -1:
                        #if fileNames.find(".csv") != -1:#push csv
                            outPutName = os.path.splitext(fileNames)[0] #this is simply a variable to name the csv output file
                            if outPutName.find("Users") != -1:
                                ViralType = outPutName.split("/")[7]
                            else:
                                ViralType = outPutName.split("/")[3]
                            print(ViralType)
                            print(outPutName+'.csv')#to follow along in the terminal


                            with open(outPutName +'.csv', 'rb') as f:
                                reader = csv.reader(f)
                                your_list = list(reader)

                            for i in  your_list:
                                if i[1].find("Accession") == -1:
                                    i.insert( 0 , ViralType )
                                    writer.writerow(i)

    return[]



protIDret()
