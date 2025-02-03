#!/usr/bin/python3

import pandas as pd

List = ["PacBio/Duplicated.tsv.gz", "PacBio/clean.tsv.gz"]

def TB_clean(File):
    TB = pd.read_csv(File, sep = '\t')
    # keep the IGH only
    TB = TB[TB.locus == 'IGH']
    # na remove
    TB = TB[~TB.v_sequence_alignment_aa.isna()]
    # remove the aa seq which has stop codon
    TB = TB[[ "*" not in i for i in TB.v_sequence_alignment_aa]]
    # remove the aa which shorter than 60
    TB = TB[[ len(i)>=60 for i in TB.v_sequence_alignment_aa]]
    # remove the seq which missed more than 4 aa on the head
    TB = TB[TB.v_germline_start<=30]
    # remove the alignment len different from the germline:
    TB= TB[[len(seq)!=len(Seq) for seq,Seq in zip(TB.v_sequence_alignment_aa, TB.v_germline_alignment_aa)]]

    TB.to_csv(File.replace('.tsv', "_IGH.tsv"), sep = '\t', index = None)


for File in List:
    TB_clean(File)
