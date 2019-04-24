import os
os.chdir('/Users/whitneyware/db/2')

infile = 'SILVA_132_SSURef_Nr99_taxonomy.rdp.outfile.txt'
outfile = 'SILVA_132_SSURef_Nr99_taxonomy.rdp.outfile_7_level_correctLabels_noSpaces.txt'

with open(infile, 'r') as f1:
    with open(outfile, 'w') as f2:
        for line in f1.readlines():
            new_line = line.replace('D_0__', 'k__')\
                .replace('D_1__', 'p__')\
                .replace('D_2__', 'c__')\
                .replace('D_3__', 'o__')\
                .replace('D_4__', 'f__')\
                .replace('D_5__', 'g__')\
                .replace('D_6__', 's__')
            new_line = new_line.split(';D_7__')
            sevenLevels = new_line[0]
            fields = sevenLevels.split('\t')
            tax = fields[1].replace(' ', '_')
            print(fields[0] + '\t' + tax, file=f2)

f1.close()
f2.close()
