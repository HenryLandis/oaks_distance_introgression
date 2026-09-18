import pandas as pd

df = pd.read_csv("tempsites.vcf", sep = '\t', header = None, on_bad_lines = 'skip')
df = df.dropna(thresh = 143, axis = 0)
df.to_csv("tempsites.vcf", sep = '\t', header = False, index = False)