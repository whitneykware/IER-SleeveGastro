infile = 'meta.txt'
outfile = 'tumorVolumeNA.txt'

sampleID = []
volume = []
d = {}

with open(infile, 'r') as file:
    with open(outfile, 'w') as out:
        for line in file.readlines():
            line = line.rstrip()
            col = line.split('\t')
            sampleID.append(col[0])
            volume.append(col[6])
        sampleID.remove('SampleID')
        volume.remove('TumorVolume')
        d = dict(zip(sampleID, volume))

        for id, vol in d.items():
            if vol == 'NA':
                print(id + '\t '+ vol, file=out)


file.close()
out.close()
