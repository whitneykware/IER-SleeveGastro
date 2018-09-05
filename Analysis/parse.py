import os
os.chdir('/Users/whitneyware/IER_Sleeve')

IER_file = 'SleeveGastro_genus_model.txt'
outfile = 'SleeveGastro_genus_Model_noMDS.txt'

with open(IER_file, 'r') as f1:
    with open(outfile, 'w') as out:
        for line in f1.readlines():
            line = line.strip()
            if not line.startswith('MDS'):
                print(line, file = out)