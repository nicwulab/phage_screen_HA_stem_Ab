import pandas as pd
from collections import Counter


TB = pd.read_csv('Result/Top100ForUpload.tsv', sep='\t')
TB2 = pd.read_csv("Result/duplicates_IGH_AA.txt", sep='\t', header=None, nrows = 100)


TB2.tail()

TB_Cmatrix = pd.DataFrame()
for i in range(100):
    Reads = TB2[1].iloc[i].split(', ')
    ID = Reads[0]
    tb_count = pd.DataFrame(Counter([i.split('-')[0] for i in Reads]), index= [0])
    tb_count['sequence_id'] = ID
    TB_Cmatrix = pd.concat([TB_Cmatrix, tb_count], ignore_index=True)

TB_Cmatrix = TB_Cmatrix.fillna(0)

TB.sequence_id.isin(TB_Cmatrix.sequence_id)
TB_result = pd.merge(TB_Cmatrix, TB, on='sequence_id', how='left')

TB_result.to_csv('Result/Top100ForUpload.tsv', sep='\t', index=False)
