#!/usr/bin/python3
import random

import numpy as np
import pandas as pd
#import logomaker
#import matplotlib.pyplot as plt
from collections import Counter
#import seaborn as sns
from sklearn.decomposition import PCA
#from umap import UMAP

#import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score  # or any other metric you wish to use
# train based on the family
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.model_selection import train_test_split


def TB_clean(TB, TB_KABAT):
    TB = TB.iloc[:,2:]
    TB_KABAT = TB_KABAT.iloc[:,2:]

    TB = TB[TB_KABAT.H1!='-']
    TB_KABAT = TB_KABAT[TB_KABAT.H1!='-']
    TB = TB.reindex(range(len(TB)))
    TB_KABAT = TB_KABAT.reindex(range(len(TB_KABAT)))

    TB = TB[TB_KABAT.H113!='-']
    TB_KABAT = TB_KABAT[TB_KABAT.H113!='-']
    TB = TB.reindex(range(len(TB)))
    TB_KABAT = TB_KABAT.reindex(range(len(TB_KABAT)))
    
    return TB, TB_KABAT

TB_KABAT = pd.read_csv('Result/Du_AA_KABAT.csv', index_col=0)
TB = pd.read_csv('Result/Du_AA_counts.csv', index_col=0)
TB['Name'] = 'ID'

# clean
TB, TB_KABAT = TB_clean(TB, TB_KABAT)

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
#TB['P23'] = TB[['P2_cm', 'P2_wt', 'P3_cm', 'P3_wt']].sum(axis=1)
TB['P23'] = TB[['P2_wt', 'P3_wt']].sum(axis=1)
#TB = TB[TB.v_family=='IGHV3-23']
TB1 = TB[TB.P23>=5]
TB2 = TB[TB.P23<5].iloc[random.sample(range(68684), k=int(len(TB)* 0.0264), ),:]


TB2 =pd.concat([TB1, TB2])
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

#Y_set = np.array(np.log10(np.array(Y_set) + 1).astype(int).tolist())
Y_set = np.array(Y_set)
#Y_set[Y_set>1] =1

X = np.array(X_set)  # 100 samples, 10 features
y = np.array([[i] for i in Y_set])  # 100 target values

# Splitting the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
X_train = np.array(X_train)

# Define the model
model = keras.Sequential([
    layers.Dense((64), activation='relu', input_shape=(X_train.shape[1],)),
    layers.Dense(64, activation='relu'),
    #layers.Dense(32, activation='relu'),
    #layers.Dense(128, activation='relu'),
    #layers.Dense(128, activation='relu'),
    layers.Dense(1)  # Output layer with one neuron for regression
])

# Compile the model
model.compile(optimizer='adam', loss='mean_squared_error')
'''
# other loss function
#from tensorflow.keras.losses import MeanSquaredLogarithmicError
from tensorflow.keras.losses import logcosh
model.compile(optimizer='adam', loss=logcosh)

# change learning rate:
from tensorflow.keras.optimizers import Adam
adam = Adam(learning_rate=0.001)
model.compile(optimizer=adam, loss='mean_squared_error')

'''
# Training the model
model.fit(X_train, y_train, epochs=100, batch_size=64, validation_split=0.2)

# Evaluating the model
loss = model.evaluate(X_test, y_test)
print(f"Test Loss: {loss}")

# Making predictions
predictions = model.predict(X_test)


def Predict_result(predictions):
    for AD in range(10):
        Test = pd.DataFrame([[i[0] for i in y_test], [i[0] for i in (AD * .1+predictions).astype(int)] ]).T

        print(AD)
        print("Total Accuracy:", 100 * len(Test[Test[0]==Test[1]])/len(Test), "%: ", len(Test))
        tmp_n = Test[Test[0]>0]
        try:
            FN =100*len(tmp_n[tmp_n[1]==0])/len(tmp_n)
        except:
            FN = 0
        print("false negative:", FN, "%: ", len(tmp_n))

        tmp_p = Test[Test[1]>0]
        try:
            FP = 100*len(tmp_p[tmp_p[0]!=0])/len(tmp_p)
        except:
            FP = 0
        print("false Positive:", FP, "% :", len(tmp_p))


def Predict_result_Yall(predictions, T = 1):
    Pre = predictions.copy()
    Pre[Pre > T ] = T
    Pre[Pre < T ] = 0
    Y_test = np.array([i[0] for i in y_test])
    Y_test[Y_test > T] = T
    Y_test[Y_test < T] = 0
    for AD in range(10):
        Test = pd.DataFrame([Y_test, [i[0] for i in (AD * .1+Pre).astype(int)] ]).T

        print(AD)
        print("Total Accuracy:", 100 * len(Test[Test[0]==Test[1]])/len(Test), "%: ", len(Test))
        tmp_n = Test[Test[0]>0]
        try:
            FN =100*len(tmp_n[tmp_n[1]==0])/len(tmp_n)
        except:
            FN = 0
        print("false negative:", FN, "%: ", len(tmp_n))

        tmp_p = Test[Test[1]>0]
        try:
            FP = 100*len(tmp_p[tmp_p[0]!=0])/len(tmp_p)
        except:
            FP = 0
        print("false Positive:", FP, "% :", len(tmp_p))




Y_test = np.array([i[0] for i in y_test])
Y_test[Y_test > T] = T
Y_test[Y_test < T] = 0
for AD in range(-20,10):
    Pre = predictions.copy()
    Pre += AD
    Pre[Pre > T ] = T
    Pre[Pre < T ] = 0
    Test = pd.DataFrame([Y_test, [i[0] for i in (Pre).astype(int)] ]).T

    print(AD)
    print("Total Accuracy:", 100 * len(Test[Test[0]==Test[1]])/len(Test), "%: ", len(Test))
    tmp_n = Test[Test[0]>0]
    try:
        FN =100*len(tmp_n[tmp_n[0]==0])/len(tmp_n)
    except:
        FN = 0
    print("false negative:", FN, "%: ", len(tmp_n))

    tmp_p = Test[Test[1]>0]
    try:
        FP = 100*len(tmp_p[tmp_p[1]!=0])/len(tmp_p)
    except:
        FP = 0
    print("false Positive:", FP, "% :", len(tmp_p))





























# TF
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.model_selection import train_test_split
import numpy as np


pca = PCA(n_components=2)
principalComponents = pca.fit_transform(TB_KABAT)
Com_TB = pd.DataFrame(pca.components_).T
Com_TB = pd.concat([Com_TB, pd.DataFrame(TB_KABAT.columns, columns=['Kabat'])], axis=1)
List = abs(Com_TB.iloc[:,:2].to_numpy().flatten())
List.sort()
List = List[::-1]
Coms = []
N_features = 40
N = 0
while len(Coms) <= N_features:
    Coms += Com_TB.Kabat[abs(Com_TB[0])==List[N]].tolist()
    Coms += Com_TB.Kabat[abs(Com_TB[1])==List[N]].tolist()
    N += 1
    Coms = list(np.unique(Coms))
Coms = Coms[:N_features]
Coms.sort()

X = TB_KABAT[Coms].copy()
#X = TB_KABAT.copy()
X = X.to_numpy().tolist()
y = TB.P3_wt.to_numpy().copy()

X_set = []
Y_set = []
for i in range(len(y)):
    X_set += [X[i]] * int(y[i]/4 +1)
    Y_set += [y[i]] * int(y[i]/4 +1)

Y_set = np.log10(np.array(Y_set) + 1).astype(int).tolist()


X = np.array(X_set)  # 100 samples, 10 features
y = np.array([[i] for i in Y_set])  # 100 target values

# Splitting the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
X_train = np.array(X_train)
# Define the model
model = keras.Sequential([
    layers.Dense(64, activation='relu', input_shape=(X_train.shape[1],)),
    layers.Dense(64, activation='relu'),
    layers.Dense(64, activation='relu'),
    layers.Dense(64, activation='relu'),
    layers.Dense(64, activation='relu'),
    layers.Dense(1)  # Output layer with one neuron for regression
])

# Compile the model
model.compile(optimizer='adam', loss='mean_squared_error')

# Training the model
model.fit(X_train, y_train, epochs=100, batch_size=32, validation_split=0.2)

# Evaluating the model
loss = model.evaluate(X_test, y_test)
print(f"Test Loss: {loss}")

# Making predictions
predictions = model.predict(X_test)

Test = pd.DataFrame([[i[0] for i in y_test], [i[0] for i in (.5+predictions).astype(int)] ]).T
len(Test[Test[0]==Test[1]])/len(Test)


model.save('my_model.h5')

with open("my_model.col", 'w') as F:
    F.write(".5 "+' '.join(Coms))

'''
from tensorflow import keras
model = keras.models.load_model('my_model.h5')
Coms = open("Result/my_model.col", 'r').read().split(' ')[1:]
model.predict(TB_KABAT[Coms].iloc[:10,:].to_numpy().tolist())
'''


# pipeline for all
X = TB_KABAT.copy()
X = X.to_numpy().tolist()
y = TB.P3_wt.to_numpy().copy()

X_set = []
Y_set = []
for i in range(len(y)):
    X_set += [X[i]] * int(np.log(y[i] + 1))
    Y_set += [y[i]] * int(np.log(y[i] + 1))

Y_set = np.log10(np.array(Y_set) + 1).astype(int).tolist()

X_train, X_test, y_train, y_test = train_test_split(X_set, Y_set, test_size=0.5, random_state=42)
dtrain = xgb.DMatrix(X_train, label=y_train)
dtest = xgb.DMatrix(X_test, label=y_test)

param = {
    'max_depth': 300,  # the maximum depth of each tree
    'eta': 0.9,      # the training step for each iteration
    'objective': 'reg:squarederror',  # regression
}
num_round = 100  # the number of training iterations
bst = xgb.train(param, dtrain, num_round)
# Make prediction
preds = bst.predict(dtest)

np.array(y[:100]).astype(int)
bst.predict(xgb.DMatrix(X[:100])).astype(int)

np.array(y[-100:]).astype(int)
bst.predict(xgb.DMatrix(X[-100:])).astype(int)

# Assuming a binary classification; adjust threshold as needed
predictions = [round(value) for value in preds]

# Evaluate accuracy
accuracy = accuracy_score(y_test, predictions)
print(f"Accuracy: {accuracy * 100.0}%")
bst.save_model('xgb_model.json')





'''
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(TB_KABAT)
Com_TB = pd.DataFrame(pca.components_).T
Com_TB = pd.concat([Com_TB, pd.DataFrame(TB_KABAT.columns, columns=['Kabat'])], axis=1)
List = abs(Com_TB.iloc[:,:2].to_numpy().flatten())
List.sort()
List = List[::-1]
Coms = []
N_features = 70
N = 0
while len(Coms) <= N_features:
    Coms += Com_TB.Kabat[abs(Com_TB[0])==List[N]].tolist()
    Coms += Com_TB.Kabat[abs(Com_TB[1])==List[N]].tolist()
    N += 1
    Coms = list(np.unique(Coms))
Coms = Coms[:N_features]



features = np.array(TB_KABAT[Coms])

umap_2d = UMAP(n_components=2, init='random', random_state=0)
proj_2d = umap_2d.fit_transform(features)


TB_U = pd.DataFrame(proj_2d, columns= ['UMAP1', "UMAP2"])
TB_U = pd.concat([TB_U, TB[['Counts', "v_family"]]], axis=1)

fg = sns.FacetGrid(TB_U, col='v_family', col_wrap=8)
fg.map(sns.scatterplot, 'UMAP1','UMAP2')
plt.show()


#umap_3d = UMAP(n_components=3, init='random', random_state=0)

#proj_3d = umap_3d.fit_transform(features)




### xgboost





family = 'IGHV1-69'
TB_matrix = TB_KABAT.copy()#[TB.v_family=="IGHV1-69"]

X = TB_matrix.copy()
X = X.astype(int)

#y = TB.P3_wt[TB.v_family=="IGHV1-69"].to_numpy().copy()
y = TB.P3_wt.to_numpy().copy()

y = np.log(np.log10(y+1).astype(int)+2).astype(int)

Head_n = 5000
X = X.head(Head_n)
y = y[:Head_n]

# PCA reduce features

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(X)
Com_TB = pd.DataFrame(pca.components_).T
Com_TB = pd.concat([Com_TB, pd.DataFrame(TB_KABAT.columns, columns=['Kabat'])], axis=1)

N_features = 40

Com_TB2 = pd.concat([Com_TB.iloc[:,[0,2]], Com_TB.iloc[:,[1,2]]])
Com_TB2[0][Com_TB2[0].isna()] = Com_TB2[1][Com_TB2[0].isna()]
Com_TB2[0] = Com_TB2[0].abs()
Com_TB2 = Com_TB2.sort_values(0, ascending = False)
Com_TB2 = Com_TB2[~Com_TB2.Kabat.duplicated()]
Coms = Com_TB2.Kabat.head(N_features).to_list()
X = X[Coms]


X = TB_KABAT.copy()
X = X.to_numpy().tolist()
y = TB.P3_wt.to_numpy().copy()

X_set = []
Y_set = []
for i in range(len(y)):
    X_set += [X[i]] * int(y[i])
    Y_set += [y[i]] * int(y[i])

Y_set = np.log10(np.array(Y_set) + 1).astype(int).tolist()

X_train, X_test, y_train, y_test = train_test_split(X_set, Y_set, test_size=0.3, random_state=42)
dtrain = xgb.DMatrix(X_train, label=y_train)
dtest = xgb.DMatrix(X_test, label=y_test)

param = {
    'max_depth': 300,  # the maximum depth of each tree
    'eta': 0.9,      # the training step for each iteration
    'objective': 'reg:squarederror',  # regression
}
num_round = 100  # the number of training iterations
bst = xgb.train(param, dtrain, num_round)
# Make prediction
preds = bst.predict(dtest)

np.array(y[:100]).astype(int)
bst.predict(xgb.DMatrix(X[:100])).astype(int)

np.array(y[-100:]).astype(int)
bst.predict(xgb.DMatrix(X[-100:])).astype(int)

# Assuming a binary classification; adjust threshold as needed
predictions = [round(value) for value in preds]

# Evaluate accuracy
accuracy = accuracy_score(y_test, predictions)
print(f"Accuracy: {accuracy * 100.0}%")
bst.save_model('xgb_model.json')
'''