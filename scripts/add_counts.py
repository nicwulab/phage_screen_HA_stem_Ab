#!/usr/bin/python3
import pandas as pd 

with open("Result/duplicates.txt", 'r') as f:
    Du_lst = f.read().split('\n')[:-1]


Du_dic = {}
[Du_dic.update({i.replace(",", '').split('\t')[-1].split(' ')[0]:i.replace(",", '').split('\t')[-1].split(' ')[1:]}) for i in Du_lst]

df = pd.read_csv('PacBio/clean.tsv.gz', compression='gzip', header=0, sep='\t')

import multiprocessing
# Assume df is the DataFrame you are working with, and Du_dic is a dictionary
def process_id(id, df, Du_dic):
    print(id)
    TMP = df[df.sequence_id == id]
    ID_du_lst = Du_dic[id]
    Rep = len(ID_du_lst)
    if not TMP.empty:
        TMP_lst = [[ii] + i[0][1:].tolist() for i, ii in zip([TMP.to_numpy()] * Rep, ID_du_lst)]
    else:
        tmp1 = df[df.sequence_id == id + "_up"]
        tmp2 = df[df.sequence_id == id + "_dn"]
        tmp_l1 = [[ii + "_up"] + i[0][1:].tolist() for i, ii in zip([tmp1.to_numpy()] * Rep, ID_du_lst)]
        tmp_l2 = [[ii + "_dn"] + i[0][1:].tolist() for i, ii in zip([tmp2.to_numpy()] * Rep, ID_du_lst)]
        TMP_lst = tmp_l1 + tmp_l2
    return TMP_lst

# This is your list of IDs
id_list = list(Du_dic.keys())

# Function to initialize the DataFrame in each worker
def init_worker(df_to_share):
    global df
    df = df_to_share

# Convert the DataFrame to a shared object to prevent copying it for each process
df_shared = df.copy()

# Create a pool of workers and distribute the tasks
if __name__ == '__main__':  # This check is necessary on Windows
    with multiprocessing.Pool(processes=64, initializer=init_worker, initargs=(df_shared,)) as pool:
        results = pool.starmap(process_id, [(id, df_shared, Du_dic) for id in id_list])

# Flatten the list of results and combine them into one list
Result = [item for sublist in results for item in sublist]
'''
# Now Result contains all the results

Result  = []
N = 0
for id in Du_dic.keys():
    N+=1
    print(N)
    TMP = df[df.sequence_id== id]
    if len(TMP) !=0:
        TMP_lst = [[ii] + i[0][1:].tolist() for i,ii in zip([TMP.to_numpy()] * Rep, ID_du_lst)]
    else:
        ID_du_lst = Du_dic[id]
        Rep = len(ID_du_lst)
        tmp1  = df[df.sequence_id== id + "_up"]
        tmp2  = df[df.sequence_id== id + "_dn"]
        tmp_l1 = [[ii+"_up"] + i[0][1:].tolist() for i,ii in zip([tmp1.to_numpy()] * Rep, ID_du_lst)]
        tmp_l2 = [[ii+"_dn"] + i[0][1:].tolist() for i,ii in zip([tmp2.to_numpy()] * Rep, ID_du_lst)]
        TMP_lst = tmp_l1 + tmp_l2
    Result +=[]



'''