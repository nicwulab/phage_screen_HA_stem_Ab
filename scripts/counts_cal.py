#!/usr/bin/python3
import pandas as pd

# Calculating the ratio matrix 

TB_du = pd.read_csv("Result/Counts_duplciate.csv", header=None, sep = ' ')
TB_du.columns = ['ID', "Sample", 'Counts']

TB_total = pd.read_csv("Result/Counts_total.csv", header=None, sep = ' ')
TB_total.columns = ['ID', "Sample", 'Counts']

TB_du['ratio'] = 0

for i in TB_du.Sample.unique():
    TB_du.ratio[TB_du.Sample == i] = TB_du.Counts[TB_du.Sample == i]/TB_total.Counts[TB_total.Sample== i].iloc[0]

TB_du.ratio*=100
TB_du = pd.concat([TB_du,pd.DataFrame([i.split('_') for i in TB_du.Sample], columns= ['Round', 'Type'])], axis= 1)

TB_du.to_csv('Result/Ratio_matrix.csv')

with open("Result/Ratio_100_list.txt", 'w') as f:
    f.write("\n".join([str(i) for i in TB_du.sort_values('ratio', ascending= False).ID.unique()[:100]]))

# Exponential Regression
import numpy as np
import multiprocessing

def Get_ratoi(ID, Type):
    TMP = TB_du[TB_du.ID==ID]
    try:
        P0 = TMP.ratio[TMP.Round=='P0'].iloc[0]
    except:
        P0 = TB_du.ratio.min()
    CM = [P0]
    WT = [P0]

    TMP_cm = TMP[TMP.Type == Type]
    for Round in ['P2', "P3"]:
        try:
            CM += [TMP_cm.ratio[TMP_cm.Round == Round].iloc[0]]
        except:
            CM += [TB_du.ratio.min()]
    return CM

ID_lst = TB_du.ID.unique()


def process_id(ID):
    results = []
    for measurement in ['cm', 'wt']:
        y = Get_ratoi(ID, measurement)
        fit = np.polyfit(x, np.log(y), 1)
        exp_fit = np.exp(fit[0])
        print(exp_fit)
        results.append([ID, measurement, exp_fit])
    return results

# This is your list of IDs
#ID_lst = range(30)  # Replace with your actual list of IDs
x = [10, 20, 30]

# Create a pool of workers equal to the number of available CPU cores
with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
    # Map process_id across all IDs in the list
    results = pool.map(process_id, ID_lst)

# Flatten the list of results
Result = [item for sublist in results for item in sublist]

# Result now contains all the results
#print(Result)
ExpRe = pd.DataFrame(Result, columns=['ID', 'Type', 'exp'])
ExpRe = ExpRe.sort_values('exp', ascending=False)
ExpRe.to_csv('Result/ExpRe.csv')
with open("Result/ExpRe_100_list.txt", 'w') as f:
    f.write("\n".join([str(i) for i in ExpRe.ID.unique()[:100]]))


'''
Result = []
x = [10, 20, 30]
for ID in ID_lst:
    y = Get_ratoi(ID, 'cm')
    fit = np.polyfit(x, np.log(y), 1)
    print(np.exp(fit[0]))
    Result += [[ID, 'cm', np.exp(fit[0])] ]
    y = Get_ratoi(ID, 'wt')
    fit = np.polyfit(x, np.log(y), 1)
    print(np.exp(fit[0]))
    Result += [[ID, 'wt', np.exp(fit[0])] ]
       

y = [1, 2, 3]
#y = [1, 10, 20]
#view the output of the model
print(fit)


X = np.array(range(x[0],x[-1]+1))
plt.scatter(x, y)
plt.plot(X, (np.exp(fit[0])**X) * np.exp(fit[1]))
plt.show()
'''