import pandas as pd

# with open('inds_rename_vcf.txt') as f:
#     colnames = [line.strip() for line in f]
    
df = pd.read_csv("tempsites.vcf", sep = '\t')
# df.columns = colnames

df_gam_gri = df[df.columns.drop(list(df.filter(regex='turbinella')))]
df_gam_gri.to_csv("pca_075_gam_gri.vcf", sep = '\t', index = False)

df_gam_tur = df[df.columns.drop(list(df.filter(regex='grisea')))]
df_gam_tur.to_csv("pca_075_gam_tur.vcf", sep = '\t', index = False)

df_gri_tur = df[df.columns.drop(list(df.filter(regex='gambelii')))]
df_gri_tur.to_csv("pca_075_gri_tur.vcf", sep = '\t', index = False)