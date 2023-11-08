#!/usr/bin/python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','-I','--input',nargs='+', help = 
                    '''input the annotated blast results. It is a list.
                    You could input the Cleared results only. But if you want to cont the ratio, please add the duplicated annotation result, too.
                    ''')   
parser.add_argument('-o','-O','--output', help= 'Please remember that the output is a directory')

# For some other parameters
parser.add_argument('-fv','-FV','--vfamily', nargs='?')
parser.add_argument('-t','-T','--top', nargs='?', type = int)
parser.add_argument('-s','-S','--sample', nargs='?')
parser.add_argument('-p','-P','--pstage', nargs='?')

args = parser.parse_args()
INPUT = args.input
OUTPUT = args.output


VFAMILY = args.vfamily
TOP = args.top
SAMPLE = args.sample
STAGE = args.pstage

import pandas as pd
import os
import numpy as np

def Read_Inputs(INPUT):
    '''
    This function is designed to take a list of file paths (INPUT), read each file as a tab-separated table, and then concatenate all these tables into a single DataFrame. The final DataFrame is then returned.
    '''
    TB = pd.DataFrame()
    for i in INPUT:
        TMP = pd.read_table(i, sep = '\t')
        TB = pd.concat([TB,TMP])
    
    TB = TB.reset_index()
    tmp = pd.DataFrame([[i.split('_')[0], i.split('_')[1].split('-')[0]] for i in TB.sequence_id], columns= ['Stage', 'Sample'])
    return pd.concat([TB, tmp], axis=1)

def Filtering(TB, filter_lst, Filter_dic):
    for i in filter_lst:
        if i not in ['TOP']:
            TB = TB[TB[i] ==  Filter_dic[i]]
    return TB

def Quality_fl(TB_f, L = 60, G_start = 10):
    TB_f = TB_f[TB_f.v_germline_start <= G_start]
    TB_f = TB_f[[len(i) >= L for i in TB_f.v_germline_alignment_aa]]
    return TB_f
    
Log = []

# 
# Collect the Filtering
Filter_dic = {"v_family": VFAMILY,
              "TOP": TOP,
              "Stage": STAGE,
              "Sample": SAMPLE}

filter_lst = [i for i in Filter_dic if Filter_dic[i] != None]

OUTPUT += "/" + "_".join([i + "_" + str(Filter_dic[i]) for i in Filter_dic]) + "/" 

if not os.path.exists(OUTPUT):
    os.makedirs(OUTPUT)

# read the table
TB = Read_Inputs(INPUT)

# filtering by the parameters
if len(filter_lst) > 0:
    TB_f = Filtering(TB, filter_lst, Filter_dic)

Log += ["Total number of locus " + TB_f.locus.iloc[0] +": " + str(len(TB[TB.locus == TB_f.locus.iloc[0]]))]
Log += ["Number of Seq after filtering by parameter: " + str(len(TB_f))]

# select appropriate result: alignment starts at first 10 bp of the germline
TB_f = Quality_fl(TB_f)
Log += ["Number of Seq after filtering by quality: " + str(len(TB_f))]
TB_f['Reads'] = [i.replace("_dn", '').replace("_up", '') for i in TB_f.sequence_id]


if i == 'TOP':
    Counts = pd.read_table("Result/ID_reads.txt", sep = ' ', header=None)
    Counts.columns = ['ID', 'Reads']
    Counts = Counts[Counts.Reads.isin(TB_f.Reads)]

    C_TB = Counts.ID.value_counts()
    C_np = Counts.to_numpy()
    Count_lst = [[C_np[i][0], C_np[i][1].split('_')[0],  C_np[i][1].split('_')[1].split('-')[0], C_np[i][1], C_TB[C_np[i][0]]] for i in range(len(C_np)) ]
    Counts = pd.DataFrame(Count_lst, columns= ['ID', 'Stage', 'Sample', 'Reads', 'Count'])



# save the results
TB_f.to_csv(OUTPUT + "Anno.tsv", sep = '\t')

# extract the sequences
TB_f1 = TB_f[TB_f.v_germline_start == 1]
TB_f2 = TB_f[TB_f.v_germline_start != 1]

# get the Max length and germline sequence
Max_l = np.max([len(seq) for seq in TB_f1.v_sequence_alignment_aa])
Seq_g = TB_f1.v_germline_alignment_aa[[len(seq) for seq in TB_f1.v_sequence_alignment_aa] == Max_l].iloc[0]


Seq1 = [[TB_f1.sequence_id.iloc[i] ,TB_f1.v_germline_alignment_aa.iloc[i]] for i in range(len(TB_f1))]
# Fill the head 
Seq2 = [[TB_f2.sequence_id.iloc[i] ,"*" * Seq_g.find(TB_f2.v_germline_alignment_aa.iloc[i][:10]) + TB_f2.v_germline_alignment_aa.iloc[i]] for i in range(len(TB_f2))]
Seq = Seq1 + Seq2

# Fill the tail
Seq = [[i[0], i[1] + ((Max_l - len(i[1])) * "*")] for i in Seq]
Seq_str = '\n'.join([">" +i[0] + "\n" + i[1] for i in Seq])

with open(OUTPUT + "aligned.fa", 'w') as f:
    f.write(Seq_str)

with open(OUTPUT + "log.txt", 'w') as f:
    f.write("\n".join(Log))

CMD = "weblogo --format PNG < " + OUTPUT + "aligned.fa  > " + OUTPUT + "aligned.png    -s large"
os.system(CMD)
