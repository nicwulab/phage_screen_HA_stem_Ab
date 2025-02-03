#!/usr/bin/python3

import pandas as pd
from abnumber import Chain
import xgboost as xgb
import numpy as np

#Fasta = "/mnt/UIUC/HA_Abs/Ab_epitope/result_old/Sarah_stem_antibodies.fasta"
Fasta = "test.fa"
bst = xgb.Booster()
bst.load_model("/mnt/UIUC/miniHA_Abs_Phage_Screen/Result/xgb_model.json")

def evaluate(lst):
    # Convert list to DMatrix for XGBoost
    dmat = xgb.DMatrix(np.array(lst).reshape(1, -1))
    # Get the prediction (you might need to adjust this based on your model)
    pred = bst.predict(dmat)
    return pred

with open(Fasta, 'r') as F:
    Seq = F.read()
#Chain(Seq, scheme='kabat')

TB = pd.read_csv("/mnt/UIUC/miniHA_Abs_Phage_Screen/Result/Du_AA_KABAT.csv", index_col= 0)
KABAT = pd.DataFrame([dict(Chain(i, scheme='kabat')) for i in  Seq.split('\n')[:-2] if ">" not in i])

KABAT.columns = [i for i in KABAT.columns.astype(str)]
KABAT = KABAT[TB.columns]
AA_dict = {
    "R": 1,
    "H": 2,
    "K": 3,
    "D": 4,
    "E": 5,
    "S": 11,
    "T": 12,
    "N": 13,
    "Q": 14,
    "C": 21,
    "G": 22,
    "P": 23,
    "A": 31,
    "V": 32,
    "I": 33,
    "L": 34,
    "M": 35,
    "F": 36,
    "Y": 37,
    "W": 38,
    "-": 0
}

KABAT_TB = KABAT.fillna('-')
KABAT_TB[KABAT_TB=='X'] = '-'
for aa in AA_dict.keys():
    KABAT_TB[KABAT_TB==aa] = AA_dict[aa] 

KABAT_TB = KABAT_TB.astype(int)