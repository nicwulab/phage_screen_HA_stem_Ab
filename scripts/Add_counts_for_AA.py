#!/usr/bin/python3
import pandas as pd
from collections import Counter
from abnumber import Chain
import multiprocessing

List_lof = "Result/duplicates_IGH_AA.txt"
TABLE = "PacBio/Duplicated_IGH.tsv.gz"
TABLE2 = "PacBio/clean_IGH.tsv.gz"

with open(List_lof, 'r') as F:
    Lst = F.readlines()

Count_lst = [[i.split('\t')[0], i.split('\t')[1].split(', ')[0]] for i in Lst]
TB = pd.DataFrame(Count_lst, columns=['Counts', 'ID'])
tmp = pd.DataFrame([Counter([i.split("-")[0] for i in i.split('\t')[-1].split(', ')]) for i in Lst])
tmp2 = pd.DataFrame([Counter([i.split("_")[0] for i in i.split('\t')[-1].split(', ')]) for i in Lst])
TB = pd.concat([TB, tmp, tmp2], axis=1)

TB_du = pd.read_csv(TABLE, sep = '\t')
TB_cl = pd.read_csv(TABLE2, sep = '\t')

tmp = TB_du[['sequence_id', 'v_family', "sequence_alignment_aa", "germline_alignment_aa"]][TB_du.sequence_id.isin(TB.ID)]
tmp2 = TB_cl[['sequence_id', 'v_family', "sequence_alignment_aa", "germline_alignment_aa"]][TB_cl.sequence_id.isin(TB.ID)]

TMP = pd.concat([tmp, tmp2])
TMP = TMP.reset_index()
TMP = TMP.iloc[:,1:]
TMP.columns = ['ID', 'v_family', "sequence_alignment_aa", "germline_alignment_aa"]

TB = pd.merge(TB, TMP, on='ID')
TB.Counts = TB.Counts.astype(int)

TB = TB.sort_values('Counts', ascending= False)
TB = TB[[len(seq)==len(Seq) for seq,Seq in zip(TB.sequence_alignment_aa, TB.germline_alignment_aa)]]
TB = TB.reset_index()
TB.to_csv('Result/Du_AA_counts.csv')

# Assuming TB_du is a DataFrame and is available globally

def process_sequence(index, sequence_alignment_aa, germline_alignment_aa):
    seq = sequence_alignment_aa
    Seq = germline_alignment_aa
    try:
        Seq_dict = dict(Chain(Seq, scheme='kabat'))
        for ii in range(len(Seq_dict)):
            Seq_dict[list(Seq_dict.keys())[ii]] = seq[ii]
    except:
        Seq_dict = {}
    return index, Seq_dict

def main(TB):
    pool = multiprocessing.Pool(processes=60)  # creates a pool of process, number of processes is equal to number of CPUs

    results = pool.starmap(process_sequence, [(TB.ID.iloc[i], TB.sequence_alignment_aa.iloc[i], TB.germline_alignment_aa.iloc[i]) for i in range(len(TB))])

    pool.close()  # close the pool to prevent any more tasks from being submitted to the pool
    pool.join()  # wait for the worker processes to exit

    # Construct the final Result dictionary from the results
    Result = {index: seq_dict for index, seq_dict in results}

    # Now Result has the processed data
    return  Result

def Result2KabatTB(Result):
    TB_KABAT = pd.DataFrame([Result[i] for i in Result.keys()])
    TB_KABAT = TB_KABAT.fillna('-')
    TB_KABAT.index = [i for i in Result.keys()]
    return TB_KABAT

'''
if __name__ == '__main__':
    Result = main(TB)
    TB_KABAT = Result2KabatTB(Result)

    TB_KABAT.to_csv('Result/Du_AA_KABAT.csv')
    Seq = "\n".join([">"+ID+ "\n" + "".join(i) for i,ID in zip(TB_KABAT.to_numpy(), TB.ID)])

    with open('Result/Dupli_aa_align.fa', 'w') as F:
        F.write(Seq)
    with open('Result/Dupli_aa_align.X.list', 'w') as F:
        F.write(' '.join([str(i) for i in TB_KABAT.columns.to_numpy()]))
'''

# clean the unreliable seqs
T_list = [len(TB_cl.sequence_alignment_aa[i]) == len(TB_cl.germline_alignment_aa[i])  for i in range(len(TB_cl)) ]
TB = TB_cl[T_list]
TB['ID'] = TB.sequence_id

if __name__ == '__main__':
    Result = main(TB)
    TB_KABAT = Result2KabatTB(Result)


    TB_KABAT.to_csv('Result/Cl_AA_KABAT.csv')
    '''
    Seq = "\n".join([">"+ID+ "\n" + "".join(i) for i,ID in zip(TB_KABAT.to_numpy(), TB.ID)])

    with open('Result/Cl_aa_align.fa', 'w') as F:
        F.write(Seq)
    with open('Result/Cl_aa_align.X.list', 'w') as F:
        F.write(' '.join([str(i) for i in TB_KABAT.columns.to_numpy()]))        
        
    ''' 
        
        
        
