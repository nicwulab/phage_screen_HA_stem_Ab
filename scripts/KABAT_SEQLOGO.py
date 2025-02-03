#!/usr/bin/python3
import numpy as np
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
from collections import Counter

TB_KABAT = pd.read_csv('Result/Du_AA_KABAT.csv', index_col=0)
TB = pd.read_csv('Result/Du_AA_counts.csv', index_col=0)


Family = 'IGHV3-20'
TB_Plot = TB_KABAT[TB.v_family == Family]#.head(50).copy()
TB_Plot[TB_Plot=="-"] = " "

crp_df = pd.DataFrame([Counter(TB_Plot.iloc[:,i]) for i in range(len(TB_Plot.columns))]).fillna(0)/len(TB_Plot)

# displays logos inline within the notebook;
# remove if using a python interpreter instead

# logomaker import
# create Logo object
fig = plt.plot()
fig, ax = plt.subplots(1,1, figsize=[21,2])
crp_logo = logomaker.Logo(crp_df,
                          shade_below=.5,
                          fade_below=.5,
                          color_scheme="weblogo_protein",  ax=ax)

# style using Logo methods
crp_logo.style_spines(visible=False)
crp_logo.style_spines(spines=['left', 'bottom'], visible=True)
crp_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
crp_logo.ax.set_ylabel("Ratio", labelpad=-1)
crp_logo.ax.xaxis.set_ticks_position('none')
crp_logo.ax.xaxis.set_tick_params(pad=-1)
crp_logo.ax.set_xticks(range(len(TB_KABAT.columns)), TB_KABAT.columns)
plt.subplots_adjust(bottom=0.4)
#plt.show()
plt.savefig('Picture/SeqLogo/Kabat_' +Family+ '_duplicated.png', dpi=300)
