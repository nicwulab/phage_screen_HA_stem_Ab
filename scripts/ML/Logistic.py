#!/usr/bin/python3
import numpy as np
import pandas as pd
#import logomaker
#import matplotlib.pyplot as plt
from collections import Counter
#import seaborn as sns
from sklearn.decomposition import PCA
#from umap import UMAP
import random
#import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score  # or any other metric you wish to use
from sklearn.linear_model import LogisticRegression


TB_KABAT = pd.read_csv('Result/Du_AA_KABAT.csv', index_col=0)
TB = pd.read_csv('Result/Du_AA_counts.csv', index_col=0)
TB['Name'] = 'ID'


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

for aa in AA_dict.keys():
    TB_KABAT[TB_KABAT==aa] = AA_dict[aa] 

TB = TB.fillna(0)
TB['P23'] = TB[['P2_cm', 'P2_wt', 'P3_cm', 'P3_wt']].sum(axis=1)

TB1 = TB[TB.P23>=5]
TB2 = TB[TB.P23<5].iloc[random.sample(range(68684), k=int(len(TB)* 0.0264), ),:]


TB2 = TB1 = pd.concat([TB1, TB2])
TB_KABAT2 = TB_KABAT.iloc[TB2.index.to_list(), :]

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(TB_KABAT2)
Com_TB = pd.DataFrame(pca.components_).T
Com_TB = pd.concat([Com_TB, pd.DataFrame(TB_KABAT2.columns, columns=['Kabat'])], axis=1)
List = abs(Com_TB.iloc[:,:2].to_numpy().flatten())
List.sort()
List = List[::-1]
Coms = []
N_features = 90
N = 0
while len(Coms) <= N_features:
    Coms += Com_TB.Kabat[abs(Com_TB[0])==List[N]].tolist()
    Coms += Com_TB.Kabat[abs(Com_TB[1])==List[N]].tolist()
    N += 1
    Coms = list(np.unique(Coms))
Coms = Coms[:N_features]
Coms.sort()

X = TB_KABAT2[Coms].copy()
#X = TB_KABAT.copy()
X = X.to_numpy().tolist()
y = TB2.P23.to_numpy().copy()

X_set = []
Y_set = []
for i in range(len(y)):
    X_set += [X[i]] #* int(y[i]/4 +1)
    Y_set += [y[i]] # * int(y[i]/4 +1)

Y_set = np.array(np.log10(np.array(Y_set) + 1).astype(int).tolist())
Y_set[Y_set>1] =1

X = np.array(X_set)  # 100 samples, 10 features
y = np.array([[i] for i in Y_set])  # 100 target values

# Splitting the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
X_train = np.array(X_train)

model = LogisticRegression()

# Fit the model with the data
model.fit(X, y)

# Example of using the model to predict a new sample
prediction = model.predict(X_test)

# Output the prediction
prediction




#Test = pd.DataFrame([[i[0] for i in y_test], (prediction[:,1]+.915).astype(int)]).T
Test = pd.DataFrame([[i[0] for i in y_test], prediction]).T

print("Total Accuracy:", 100 * len(Test[Test[0]==Test[1]])/len(Test), "%: ", len(Test))
tmp_n = Test[Test[0]>0]
print("false negative:", 100*len(tmp_n[tmp_n[1]==0])/len(tmp_n), "%: ", len(tmp_n))

tmp_p = Test[Test[1]>0]
print("false Positive:", 100*len(tmp_p[tmp_p[0]!=0])/len(tmp_p), "% :", len(tmp_p))




