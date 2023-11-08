import pandas as pd
import numpy as np
import seqlogo


TB = pd.read_table("Result/HV1-69/HV1-69_domain_all.tsv", sep ='\t')
TB = TB[TB.v_germline_start <= 1]
#TB2 = TB[TB.v_germline_start.isin([1,3,5])]

Seq = TB[[  "sequence_alignment_aa", 
            "v_germline_start"]].to_numpy()

Seq_1 = [Seq[i][0][:120]  for i in range(len(Seq))]
Seq_1 = [Seq_1[i] + ("*" * (120 - len(Seq_1[i])))  for i in range(len(Seq_1))]
with open('test.fa', 'w') as F: F.write("\n".join([">1\n"+i for i in Seq_1]))