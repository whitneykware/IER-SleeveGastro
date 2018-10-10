import os
os.chdir('/Users/whitneyware/IER_Sleeve')

IER_file = 'IER_genus_model_noMDS.txt'
SG_file = 'SleeveGastro_genus_model_noMDS.txt'
outfile = 'merged_pVals.txt'

IER_bug = []
SG_bug = []
IERpVal_time = []
IERpVal_tumor = []
IERpVal_group = []
IERcorr = []
IERlog10pValDir = []

SGpVal_time = []
SGpVal_tumor = []
SGpVal_group = []
SGcorr = []
SGlog10pValDir = []

with open(IER_file, 'r') as f1:
    with open(SG_file, 'r') as f2:
        for line in f1.readlines():
            if line.startswith('bugNames'):
                continue
            line = line.strip()
            col = line.split('\t')
            IER_bug.append(col[0])
            IERpVal_time.append(col[1])
            IERpVal_tumor.append(col[2])
            IERpVal_group.append(col[3])
            IERlog10pValDir.append(col[9])
            IERcorr.append(col[11])

        for line in f2.readlines():
            if line.startswith('bugNames'):
                continue
            line = line.strip()
            col = line.split('\t')
            SG_bug.append(col[0])
            SGpVal_time.append(col[1])
            SGpVal_tumor.append(col[2])
            SGpVal_group.append(col[3])
            SGlog10pValDir.append(col[6])
            SGcorr.append(col[7])
f1.close(), f2.close()

#bugs = IER_bug + SG_bug
both = set(IER_bug).intersection(SG_bug)

IER_time = dict(zip(IER_bug,IERpVal_time))
IER_tumor = dict(zip(IER_bug,IERpVal_tumor))
IER_group = dict(zip(IER_bug,IERpVal_group))
IER_log10 = dict(zip(IER_bug, IERlog10pValDir))
IER_corr = dict(zip(IER_bug, IERcorr))

SG_time = dict(zip(SG_bug,SGpVal_time))
SG_tumor = dict(zip(SG_bug,SGpVal_tumor))
SG_group = dict(zip(SG_bug,SGpVal_group))
SG_log10 = dict(zip(SG_bug, SGlog10pValDir))
SG_corr = dict(zip(SG_bug, SGcorr))


with open(outfile, 'w') as out:
    print("taxa\tIERGroupPValue\tIERTimePValue\tIERTumorPValue\tIERpValuesLog10VolumeWithDirection\t"
          "IERcorCoeffVolumeBUG\tSleeveGroupPValue\tSleeveTimePValue\tSleeveTumorPValue\t"
          "SleevePValuesLog10VolumeWithDirection\tSleeveCorCoeffVolumeBUG", file=out)
    for i in both:
        print(i +'\t' + IER_group[i]+'\t' + IER_time[i]+'\t' + IER_tumor[i]+'\t' + IER_log10[i] + '\t' +
              IER_corr[i]+ '\t' + SG_group[i]+'\t' + SG_time[i]+'\t' + SG_tumor[i]+'\t' +SG_log10[i]+ '\t' +
              SG_corr[i], file=out)
out.close()