import re
import os

#os.chdir('/Users/whitneyware/IER_Project/seqs')
seqInfile = 'combined_seqs.fna'
#allOut = 'combined_allSampleID2.txt'
finalOut = 'combined_seqSampleID.txt'

sampleID = []
with open(seqInfile, 'r') as file:
    #with open(allOut, 'w') as out:
        with open(finalOut, 'w') as outf:
            for line in file.readlines():
                line = line.strip()
                if line.startswith('>'):
                    header = line.split('_')
                    ID = header[0]
                    ID = re.sub('(>)', '', ID)
                    sampleID.append(ID)

            #for x in sampleID:
             #   print(x, file=out)
            uni = set(sampleID)
            for i in uni:
                print(i, file=outf)


file.close()
#out.close()
outf.close()