import os
os.chdir('/Users/whitneyware/IER_Sleeve/IER_Run_2_seqs')

barcodeLen = 14

with open('barcodes.fastq', 'r') as bc:
    with open("IER_run_2_corrected_barcodes.fastq", 'w') as out:
        for line in bc.readlines():
            line = line.rstrip()
            if line.startswith("@" or "+" or "'"):
                print(line, file=out)
            else:
                line = line[:barcodeLen]
                print(line, file=out)

bc.close()
out.close()
