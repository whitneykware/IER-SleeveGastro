import os

os.chdir("/Users/whitneyware/IER_Sleeve")

seq = []
taxa = []
d = {}

with open('dada2_IER_taxa_names.txt', 'r') as t:
    #with open('dada2_ier_dic.txt', 'w') as out:
        for line in t.readlines():
            line = line.rstrip()
            cols = line.split('\t')
            seq = cols[0]
            taxa = cols[1]
            d[seq] = taxa


t.close()

with open("IER_fwdOnly_DADA2_Model.txt", "r") as m:
    with open("IER_taxa_dada2_model.txt", "w") as o:
        text = m.read()

        for key in d:
            text = text.replace(key, d[key])

        o.write(text)

m.close()
o.close()
