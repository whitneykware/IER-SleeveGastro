import os

os.chdir("/Users/whitneyware/IER_Sleeve/sleeve_seqs")


with open('sleeve_tax_tab.txt', 'r') as t:
    with open('dada2_sleeve_taxa_names.txt', 'w') as out:
        for line in t.readlines():
            line = line.rstrip()
            if line.startswith('#'):
                continue
            cols = line.split('\t')
            seq = cols[0]
            kin = cols[1]
            phy = cols[2]
            cla = cols[3]
            order = cols[4]
            fam = cols[5]
            gen = cols[6]
            taxa_names = str(kin + ', ' + phy + ', ' + cla + ', ' + order + ', ' + fam + ', ' + gen)
            print(seq+'\t'+taxa_names, file=out)


t.close()
out.close()