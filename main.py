from Bio import SeqIO
from pathlib import Path
import os

pathlist = Path('./pass').glob('*.fastq')
forward = {'ACGAGAC': ['fw_1'], 'AATCAGT': ['fw_1'], 'CGAGACT': ['fw_1'], 'ATCAGTC': ['fw_1'], 'GAGACTG': ['fw_1'], 
 'TCAGTCT': ['fw_1'], 'AGACTGA': ['fw_1'], 'CAGTCTC': ['fw_1'], 'GACTGAT': ['fw_1'], 'AGTCTCG': ['fw_1'], 
 'GCTGTAC': ['fw_2'], 'AATCCGT': ['fw_2'], 'CTGTACG': ['fw_2'], 'ATCCGTA': ['fw_2'], 'TGTACGG': ['fw_2'], 
 'TCCGTAC': ['fw_2'], 'GTACGGA': ['fw_2'], 'CCGTACA': ['fw_2'], 'TACGGAT': ['fw_2'], 'CGTACAG': ['fw_2'], 
 'ATCACCA': ['fw_3'], 'ACACCTG': ['fw_3'], 'TCACCAG': ['fw_3'], 'CACCTGG': ['fw_3'], 'CACCAGG': ['fw_3'], 
 'ACCTGGT': ['fw_3'], 'ACCAGGT': ['fw_3'], 'CCTGGTG': ['fw_3'], 'CCAGGTG': ['fw_3'], 'CTGGTGA': ['fw_3'], 
 'TGGTCAA': ['fw_4'], 'TATCGTT': ['fw_4'], 'GGTCAAC': ['fw_4'], 'ATCGTTG': ['fw_4'], 'GTCAACG': ['fw_4'], 
 'TCGTTGA': ['fw_4'], 'TCAACGA': ['fw_4'], 'CGTTGAC': ['fw_4'], 'CAACGAT': ['fw_4'], 'GTTGACC': ['fw_4'], 
 'ATCGCAC': ['fw_5'], 'TTACTGT': ['fw_5'], 'TCGCACA': ['fw_5'], 'TACTGTG': ['fw_5'], 'CGCACAG': ['fw_5'], 
 'ACTGTGC': ['fw_5'], 'GCACAGT': ['fw_5'], 'CTGTGCG': ['fw_5'], 'CACAGTA': ['fw_5'], 'TGTGCGA': ['fw_5'], 
 'GTCGTGT': ['fw_6'], 'AGGCTAC': ['fw_6'], 'TCGTGTA': ['fw_6'], 'GGCTACA': ['fw_6'], 'CGTGTAG': ['fw_6'], 
 'GCTACAC': ['fw_6'], 'GTGTAGC': ['fw_6'], 'CTACACG': ['fw_6'], 'TGTAGCC': ['fw_6'], 'TACACGA': ['fw_6', 'rv_8'], 
 'AGCGGAG': ['fw_7'], 'CTAACCT': ['fw_7'], 'GCGGAGG': ['fw_7'], 'TAACCTC': ['fw_7'], 'CGGAGGT': ['fw_7'], 
 'AACCTCC': ['fw_7'], 'GGAGGTT': ['fw_7'], 'ACCTCCG': ['fw_7'], 'GAGGTTA': ['fw_7'], 'CCTCCGC': ['fw_7'], 
 'ATCCTTT': ['fw_8'], 'GAACCAA': ['fw_8'], 'TCCTTTG': ['fw_8'], 'AACCAAA': ['fw_8'], 'CCTTTGG': ['fw_8'], 
 'ACCAAAG': ['fw_8'], 'CTTTGGT': ['fw_8'], 'CCAAAGG': ['fw_8'], 'TTTGGTT': ['fw_8'], 'CAAAGGA': ['fw_8'], }

reverse =  {'TACAGCG': ['rv_1'], 'GTATGCG': ['rv_1'], 'ACAGCGC': ['rv_1'], 'TATGCGC': ['rv_1'], 'CAGCGCA': ['rv_1'], 
 'ATGCGCT': ['rv_1'], 'AGCGCAT': ['rv_1'], 'TGCGCTG': ['rv_1'], 'GCGCATA': ['rv_1'], 'GCGCTGT': ['rv_1'], 
 'ACCGGTA': ['rv_2'], 'GTACATA': ['rv_2'], 'CCGGTAT': ['rv_2'], 'TACATAC': ['rv_2'], 'CGGTATG': ['rv_2'], 
 'ACATACC': ['rv_2'], 'GGTATGT': ['rv_2'], 'CATACCG': ['rv_2'], 'GTATGTA': ['rv_2'], 'ATACCGG': ['rv_2'], 
 'AATTGTG': ['rv_3'], 'TCCGACA': ['rv_3'], 'ATTGTGT': ['rv_3'], 'CCGACAC': ['rv_3'], 'TTGTGTC': ['rv_3'], 
 'CGACACA': ['rv_3'], 'TGTGTCG': ['rv_3'], 'GACACAA': ['rv_3'], 'GTGTCGG': ['rv_3'], 'ACACAAT': ['rv_3'], 
 'TGCATAC': ['rv_4'], 'CCAGTGT': ['rv_4'], 'GCATACA': ['rv_4'], 'CAGTGTA': ['rv_4'], 'CATACAC': ['rv_4'], 
 'AGTGTAT': ['rv_4'], 'ATACACT': ['rv_4'], 'GTGTATG': ['rv_4'], 'TACACTG': ['rv_4'], 'TGTATGC': ['rv_4'], 
 'AGTCGAA': ['rv_5'], 'CCTCGTT': ['rv_5'], 'GTCGAAC': ['rv_5'], 'CTCGTTC': ['rv_5'], 'TCGAACG': ['rv_5'], 
 'TCGTTCG': ['rv_5'], 'CGAACGA': ['rv_5'], 'CGTTCGA': ['rv_5'], 'GAACGAG': ['rv_5'], 'GTTCGAC': ['rv_5'], 
 'ACCAGTG': ['rv_6'], 'TGAGTCA': ['rv_6'], 'CCAGTGA': ['rv_6'], 'GAGTCAC': ['rv_6'], 'CAGTGAC': ['rv_6'], 
 'AGTCACT': ['rv_6'], 'AGTGACT': ['rv_6'], 'GTCACTG': ['rv_6'], 'GTGACTC': ['rv_6'], 'TCACTGG': ['rv_6'], 
 'GAATACC': ['rv_7'], 'GACTTGG': ['rv_7'], 'AATACCA': ['rv_7'], 'ACTTGGT': ['rv_7'], 'ATACCAA': ['rv_7'], 
 'CTTGGTA': ['rv_7'], 'TACCAAG': ['rv_7'], 'TTGGTAT': ['rv_7'], 'ACCAAGT': ['rv_7'], 'TGGTATT': ['rv_7'], 
 'GTAGATC': ['rv_8'], 'TAGATCG': ['rv_8'], 'ACACGAT': ['rv_8'], 'AGATCGT': ['rv_8'], 'CACGATC': ['rv_8'], 
 'GATCGTG': ['rv_8'], 'ACGATCT': ['rv_8'], 'ATCGTGT': ['rv_8'], 'CGATCTA': ['rv_8'], 'TAACGTG': ['rv_9'], 
 'GCACACA': ['rv_9'], 'AACGTGT': ['rv_9'], 'CACACAC': ['rv_9'], 'ACGTGTG': ['rv_9'], 'ACACACG': ['rv_9'], 
 'CGTGTGT': ['rv_9'], 'CACACGT': ['rv_9'], 'GTGTGTG': ['rv_9'], 'ACACGTT': ['rv_9'], 'CATTATG': ['rv_10'], 
 'CACGCCA': ['rv_10'], 'ATTATGG': ['rv_10'], 'ACGCCAT': ['rv_10'], 'TTATGGC': ['rv_10'], 'CGCCATA': ['rv_10'], 
 'TATGGCG': ['rv_10'], 'GCCATAA': ['rv_10'], 'ATGGCGT': ['rv_10'], 'CCATAAT': ['rv_10'], 'CCAATAC': ['rv_11'], 
 'CAGGCGT': ['rv_11'], 'CAATACG': ['rv_11'], 'AGGCGTA': ['rv_11'], 'AATACGC': ['rv_11'], 'GGCGTAT': ['rv_11'], 
 'ATACGCC': ['rv_11'], 'GCGTATT': ['rv_11'], 'TACGCCT': ['rv_11'], 'CGTATTG': ['rv_11'], 'GATCTGC': ['rv_12'], 
 'GGATCGC': ['rv_12'], 'ATCTGCG': ['rv_12'], 'GATCGCA': ['rv_12'], 'TCTGCGA': ['rv_12'], 'ATCGCAG': ['rv_12'], 
 'CTGCGAT': ['rv_12'], 'TCGCAGA': ['rv_12'], 'TGCGATC': ['rv_12'], 'CGCAGAT': ['rv_12']}

start = 5
length = 50
size = 7

final = {}

for path in pathlist:
    for seq_record in SeqIO.parse(path, 'fastq'):
        for i in range(length):
            fw = str(seq_record.seq)[(start + i) : (start + size + i)]
            rv = str(seq_record.seq)[(-start -size - i) : (-start - i)]

            if fw in forward:
                seq = str(seq_record.seq)
                if seq not in final:
                    final[seq] = [forward[fw]]
                else:
                    final[seq].append(forward[fw])
            if rv in reverse:
                seq = str(seq_record.seq)
                if seq not in final:
                    final[seq] = [reverse[rv]]
                else:
                    final[seq].append(reverse[rv])
    print(final)
    print(len(final))
    break