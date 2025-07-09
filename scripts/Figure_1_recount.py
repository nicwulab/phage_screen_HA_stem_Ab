'''
This script is for Figure 1 in the manuscript.
Step 1:
The script would integrate the productive sequences from the PacBio data,
Step 2:
Then, only full productive sequences with both heavy and light chains are kept.
Step 3:
The script would pivot the table to get the heavy and light chains, and extract the screening information from the barcode.
Finally, it would count the number of sequences for each combination of screening, heavy chain, and light chain.

Step 4:
    Using the well paired heavy and light chains, to make a statistical table.
'''
import pandas as pd
from collections import Counter
import numpy as np

def Productive_clean(TB):
    TB = TB[TB.stop_codon == 'F']
    TB = TB[TB.productive == 'T']
    TB = TB[~TB.v_family.isna()]
    return TB

# Step 1: Data read and clean
# Read the data
TB1 = pd.read_csv("PacBio/clean.tsv.gz", sep="\t")
TB2 = pd.read_csv("PacBio/Duplicated.tsv.gz", sep="\t")

# Clean the data (remove non-productive sequences and sequences without v_family)
TB1 = Productive_clean(TB1)
TB2 = Productive_clean(TB2)

# merge into one table
columns = ['sequence_id', 'v_family']
tb_anno_du = pd.concat([TB1[columns], TB2[columns]], axis = 0)

# Step 2: Find the sequences with both heavy and light chains 
# Check the duplicated barcode (sequence_id) to find the sequences with H+L
sequence_id = tb_anno_du.sequence_id.to_list()
barcode = [i[:-3] for i in sequence_id]
tb_anno_du['barcode'] = barcode
barcode_full = tb_anno_du.barcode[tb_anno_du.barcode.duplicated()].to_list()
len(barcode_full)
tb_anno_full = tb_anno_du[tb_anno_du.barcode.isin(barcode_full)]
tb_anno_full['Chain'] = [i[-2:] for i in tb_anno_full.sequence_id]

# Step 3: Pivot the table to get the heavy and light chains
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



# Step 4: Statistical table

columns = ['sequence_id', 'v_family', 'sequence_alignment', 'sequence_alignment_aa']
tb2 = pd.concat([TB1[columns], TB2[columns]], axis = 0)
sequence_id = tb2.sequence_id.to_list()
barcode = [i[:-3] for i in sequence_id]
tb2['barcode'] = barcode
tb2 = tb2[tb2.barcode.isin(tb_anno_full.barcode)]
tb2['Chain'] = [i[-2:] for i in tb2.sequence_id]
tb2 = tb2[[i in ['up', 'dn'] for i in tb2['Chain']]]
tb2['Chain'] = "l" 
tb2['Chain'][["IGH" in i for i in tb2.v_family]] = 'h'

tb2_l = tb2[tb2.Chain == 'l']
tb2_h = tb2[tb2.Chain == 'h']

tb_m = pd.merge(tb2_l, tb2_h, on = 'barcode', suffixes = ('_l', '_h'))
tb_m['seq_full'] = tb_m.sequence_alignment_aa_l + tb_m.sequence_alignment_aa_h

# save all counts
Seq_Ct = pd.DataFrame(Counter(tb_m.seq_full), index = ['Count']).T
Seq_Ct.sort_values(by = 'Count', ascending = False, inplace = True)
Seq_Ct.to_csv("Result/Paired_aa_counts.tsv", sep = "\t", index=False)


Top_10 = [i[0] for i in  Counter(tb_m.seq_full).most_common(10)]  # Check how many sequences are duplicated
tb_mf = tb_m[tb_m.seq_full.isin(Top_10)]

tb_top10 = pd.DataFrame()
for aa in Top_10:
    bc = tb_mf.barcode[tb_mf.seq_full == aa].to_list()
    tb_bc = pd.DataFrame(Counter([i.split('-')[0] for i in bc]), index = ['Count'])
    tb_bc['seq_full'] = aa
    tb_top10 = pd.concat([tb_top10, tb_bc], axis = 0)

# fill na with 0
tb_top10.fillna(0, inplace=True)
tb_top10['Name']  = np.nan 

m_HB31 = ["EVQLLESGGGLVQPGGSLRLSCTVSGFTSSNYDMSWVRQAPGKGLEWVSGISGTGRDTYYADSVKGRFTISRDNSKNTLSLQMNSLRADDTAVYFCARDRPVLRYFDWQPYGLDVWGKGTTVTVSS" in  i and "SASGSPGQSVTISCTGTSSDVGGYNYVSWYQQYPSKAPKLMIYEVSKRPSGVPDRFSGSKSGNTASLTVSGLQAEDEADYYCNSYAGSNLYVFGTGTKVTVL" in  i  for i in tb_top10.seq_full]
m_HB34 = ["EVQLLESGGGLVQPGGSLRLSCTVSGFTSSNYDMSWVRQAPGKGLEWVSGISGTGRDTYYADSVKGRFTISRDNSKNTLSLQMNSLRADDTAVYFCARDRPVLRYFDWQPYGLDVWGKGTTVTVSS" in  i and "QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQHPGKAPKLMIYEGSQRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCCSYASSTSYVFGTGTKVTVL" in  i  for i in tb_top10.seq_full]
m_HB315 = ["EVQLLESGGGLAQPGGSLRLSCTASGFTFSNFAMNWVRQAPGKGLEWVSFINSGSGSTYYADSVKGRFTISRDNSENTLFLQMNSLRAEDTAVYYCAKDLGRAIFGLVIPEGAFDIWGQGTMVTVSS" in  i and "EIVMTQSPSSLSASVGDRVTITCRASQSINNHLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPADFATYHCQQSYGVPLTFGQGTKVEIK" in  i  for i in tb_top10.seq_full]


tb_top10.Name[m_HB31] = 'HB31'
tb_top10.Name[m_HB34] = 'HB34'
tb_top10.Name[m_HB315] = 'HB315'
tb_mfq = tb_mf[~tb_mf.seq_full.duplicated()]

# merge tb_mfq and tb_top10 based on column 'seq_full':
tb_top10 = tb_top10.reset_index()
tb_mfq = tb_mfq.reset_index()
tb_final = pd.merge(tb_top10, tb_mfq, on='seq_full', how='left')   

tb_final.drop(columns=['index_x', 'index_y', 'sequence_id_l', 'sequence_id_h', 'Chain_l', 'Chain_h', 'seq_full', 'P2_cm', 'P2_wt'    ], inplace=True)

tb_final.to_csv("Result/Top10_sequences.tsv", sep="\t", index=False)



