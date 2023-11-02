#!/usr/bin/python3
from Bio import SeqIO


with open("blast.out", 'r') as f:
    TB = f.read().split("\n")[:-1]
TB_dic = {}
[TB_dic.update({i.split('\t')[1]:i.split('\t')[-2:]}) for i in TB]


Seq1="PacBio/clean.fa"

Result = []
for seq_record in SeqIO.parse(Seq1, "fasta"):
    try:
        tmp = [int(i) for i in TB_dic[seq_record.id]]
        tmp.sort()
        Seq_A = "".join([
            ">", seq_record.id, "_up\n", str(seq_record.seq[:tmp[0]]), 
            "\n>", seq_record.id, "_dn\n", str(seq_record.seq[tmp[1]:])])
    except:
        Seq_A = "".join([">", seq_record.id, "\n", str(seq_record.seq)])
    Result += [Seq_A]

with open("PacBio/clean_split.fa", 'w') as f:
    f.write("\n".join(Result))
