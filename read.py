from Bio import SeqIO
from pathlib import Path
import os

pathlist = Path('./pass').glob('*.fastq')
os.mkdir('txt')

for path in pathlist:
    list_of_target_seqs = []
    for seq_record in SeqIO.parse(path, 'fastq'):
        list_of_target_seqs.append(str(seq_record.seq))

    with open(f'txt/{str(path)[17:-6]}.txt', 'w') as file:
        for seq in list_of_target_seqs:
            file.write(seq + '\n'+'\n')