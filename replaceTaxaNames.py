import os

os.chdir("/Users/whitneyware/IER_Sleeve/sleeve_seqs")


with open('sleeve_tax_tab.txt', 'r') as t:
    with open('dada2_sleeve_taxa_names.txt', 'r') as names: