##This python script will have the goal of opening a csv file containing all the
##Go terms for all genes which were 1. annotated and 2. had >=500 hits
##I think my strategy will be to open the file and to count the number of times
##each term is found in the list without regards to which gene it is associated
##to

import csv, os
from collections import Counter

os.chdir('D:/bioinformatics/research project')

def GOTermCount() :
    CleanTermList = []

    with open('GoTermOnly.csv', 'rb') as f:
        reader = csv.reader(f)
        TermList = list(reader)

        for row in TermList:

            TempString = row[6].strip('[]')
            TempString = TempString.replace("transcription," , "transcription")
            TempString = TempString.replace("activity," , "activity")

            for Terms in TempString.split(','):

                Terms = Terms.lstrip()
                CleanTermList.append(Terms.replace("'","").replace('"',''))

        AllTermCount = Counter(CleanTermList)

    CleanTermList = []


    with open('HighGO.csv', 'rb') as f:
        reader = csv.reader(f)
        TermList = list(reader)

        for row in TermList:

            TempString = row[6].strip('[]')
            TempString = TempString.replace("transcription," , "transcription")
            TempString = TempString.replace("activity," , "activity")

            for Terms in TempString.split(','):

                Terms = Terms.lstrip()
                CleanTermList.append(Terms.replace("'","").replace('"',''))

        TermCount = Counter(CleanTermList)

    matchingGo = [[k, AllTermCount[k],TermCount[k]] for k in TermCount if k in AllTermCount]

    with open('CountAllGO.csv','wb') as p:
        writer = csv.writer(p)
        writer.writerow(["Term", "allCount", "HighCount"])
        for i in matchingGo:
            writer.writerow(i)


    return[]

GOTermCount()
