import numpy as np
import pandas as pd
from scipy import stats

df_payson = pd.read_csv("dtrios_output/inds_payson_tree.txt", sep = '\t')
payson_p = df_payson['p-value'].to_numpy()
payson_fdr_p = stats.false_discovery_control(payson_p)
df_payson['p-value_FDR'] = payson_fdr_p.tolist()
df_payson.to_csv("dtrios_output/payson/inds_payson_tree_FDR.csv", index = False)

df_gila = pd.read_csv("dtrios_output/inds_gila_tree.txt", sep = '\t')
gila_p = df_gila['p-value'].to_numpy()
gila_fdr_p = stats.false_discovery_control(gila_p)
df_gila['p-value_FDR'] = gila_fdr_p.tolist()
df_gila.to_csv("dtrios_output/gila/inds_gila_tree_FDR.csv", index = False)

df_globe = pd.read_csv("dtrios_output/inds_globe_tree.txt", sep = '\t')
globe_p = df_globe['p-value'].to_numpy()
globe_fdr_p = stats.false_discovery_control(globe_p)
df_globe['p-value_FDR'] = globe_fdr_p.tolist()
df_globe.to_csv("dtrios_output/globe/inds_globe_tree_FDR.csv", index = False)

df_lincoln = pd.read_csv("dtrios_output/inds_lincoln_tree.txt", sep = '\t')
lincoln_p = df_lincoln['p-value'].to_numpy()
lincoln_fdr_p = stats.false_discovery_control(lincoln_p)
df_lincoln['p-value_FDR'] = lincoln_fdr_p.tolist()
df_lincoln.to_csv("dtrios_output/lincoln/inds_lincoln_tree_FDR.csv", index = False)

df_goldenvalley = pd.read_csv("dtrios_output/inds_goldenvalley_tree.txt", sep = '\t')
goldenvalley_p = df_goldenvalley['p-value'].to_numpy()
goldenvalley_fdr_p = stats.false_discovery_control(goldenvalley_p)
df_goldenvalley['p-value_FDR'] = goldenvalley_fdr_p.tolist()
df_goldenvalley.to_csv("dtrios_output/goldenvalley/inds_goldenvalley_FDR.csv", index = False)

df_santafe = pd.read_csv("dtrios_output/inds_santafe_tree.txt", sep = '\t')
santafe_p = df_santafe['p-value'].to_numpy()
santafe_fdr_p = stats.false_discovery_control(santafe_p)
df_santafe['p-value_FDR'] = santafe_fdr_p.tolist()
df_santafe.to_csv("dtrios_output/santafe/inds_santafe_FDR.csv", index = False)

df_sanjuan = pd.read_csv("dtrios_output/inds_sanjuan_tree.txt", sep = '\t')
sanjuan_p = df_sanjuan['p-value'].to_numpy()
sanjuan_fdr_p = stats.false_discovery_control(sanjuan_p)
df_sanjuan['p-value_FDR'] = sanjuan_fdr_p.tolist()
df_sanjuan.to_csv("dtrios_output/sanjuan/inds_sanjuan_FDR.csv", index = False)