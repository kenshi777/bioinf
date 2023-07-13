from Bio.Seq import Seq

primers = {
    'fw_1' : 'ACGAGACTGATT',
    'fw_2' : 'GCTGTACGGATT',
    'fw_3' : 'ATCACCAGGTGT',
    'fw_4' : 'TGGTCAACGATA',
    'fw_5' : 'ATCGCACAGTAA',
    'fw_6' : 'GTCGTGTAGCCT',
    'fw_7' : 'AGCGGAGGTTAG',
    'fw_8' : 'ATCCTTTGGTTC',
    'rv_1' : 'TACAGCGCATAC',
    'rv_2' : 'ACCGGTATGTAC',
    'rv_3' : 'AATTGTGTCGGA',
    'rv_4' : 'TGCATACACTGG',
    'rv_5' : 'AGTCGAACGAGG',
    'rv_6' : 'ACCAGTGACTCA',
    'rv_7' : 'GAATACCAAGTC',
    'rv_8' : 'GTAGATCGTGTA',
    'rv_9' : 'TAACGTGTGTGC',
    'rv_10' : 'CATTATGGCGTG',
    'rv_11' : 'CCAATACGCCTG',
    'rv_12' : 'GATCTGCGATCC'
}

l = 12
size = 7
tmp = []

d = {}

for num, primer in primers.items():
    for j in range(l - size):
        part = primer[j:size + j]
        rc_part = str(Seq(primer).reverse_complement())[j:size + j]
        if part not in d:
            d[part] = [num]
        else:
            d[part].append(num)
        if rc_part not in d:
            d[rc_part] = [num]
        else:
            d[rc_part].append(num)

print(d)
with open('window.txt', 'w') as file:
    file.write(str(d))

