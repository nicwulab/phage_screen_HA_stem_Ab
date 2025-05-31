import pandas as pd
from collections import Counter

def Productive_clean(TB):
    TB = TB[TB.stop_codon == 'F']
    TB = TB[TB.productive == 'T']
    TB = TB[~TB.v_family.isna()]
    return TB

# Read the data
TB1 = pd.read_csv("PacBio/clean.tsv.gz", sep="\t")
TB2 = pd.read_csv("PacBio/Duplicated.tsv.gz", sep="\t")

# Clean the data
TB1 = Productive_clean(TB1)
TB2 = Productive_clean(TB2)

# merge into one table
columns = ['sequence_id', 'v_family']
tb_anno_du = pd.concat([TB1[columns], TB2[columns]], axis = 0)

# Check the duplicated barcode
sequence_id = tb_anno_du.sequence_id.to_list()
barcode = [i[:-3] for i in sequence_id]
tb_anno_du['barcode'] = barcode
barcode_full = tb_anno_du.barcode[tb_anno_du.barcode.duplicated()].to_list()
len(barcode_full)
tb_anno_full = tb_anno_du[tb_anno_du.barcode.isin(barcode_full)]

tb_anno_full['Chain'] = [i[-2:] for i in tb_anno_full.sequence_id]

TB_plot = tb_anno_full[['barcode', 'Chain', 'v_family']].pivot(index = 'barcode', columns = 'Chain', values = 'v_family').reset_index()
# the parameter of aggfunc is not important, since we only have one value for each barcode and chain, so we can use any function here.
TB_plot = TB_plot[['barcode', 'dn', 'up']]
TB_plot = TB_plot[~TB_plot.dn.isna()]
TB_plot = TB_plot[~TB_plot.up.isna()]

TB_plot['IGH'] = pd.NA
TB_plot['IGL'] = pd.NA

mask_dn = ["IGHV" in i for i in TB_plot.dn]
mask_up = ["IGHV" in i for i in TB_plot.up]

mask_dnr = ["IGHV" not in i for i in TB_plot.dn]
mask_upr = ["IGHV" not in i for i in TB_plot.up]

TB_plot.IGH[mask_dn] = TB_plot.dn[mask_dn]
TB_plot.IGH[mask_up] = TB_plot.up[mask_up]

TB_plot.IGL[mask_dnr] = TB_plot.dn[mask_dnr]
TB_plot.IGL[mask_upr] = TB_plot.up[mask_upr]

# we need to extract the screening information from the barcode
TB_screen = pd.DataFrame([i.split('-')[0] for i in TB_plot.barcode], columns = ['P']).reset_index(drop = True)

TB_plot = TB_plot.reset_index(drop = True)
TB_plot = pd.concat([TB_plot, TB_screen], axis = 1)

pd.DataFrame(Counter(TB_plot.P + ":" + TB_plot.IGH + ":" + TB_plot.IGL), index = ['Count']).T.to_csv("Result/IGH_IGL.tsv", sep = "\t")
