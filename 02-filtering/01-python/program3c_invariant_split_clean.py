import pandas as pd

df = pd.read_csv(
    "all_sites_pieces/all_sites_01.vcf", sep = '\t', header = None, on_bad_lines = 'skip')
df = df.dropna(thresh = 143, axis = 0)
df.to_csv("all_sites_pieces/all_sites_01.vcf", sep = '\t', header = False, index = False)

df = pd.read_csv(
    "all_sites_pieces/all_sites_02.vcf", sep = '\t', header = None, on_bad_lines = 'skip')
df = df.dropna(thresh = 143, axis = 0)
df.to_csv("all_sites_pieces/all_sites_02.vcf", sep = '\t', header = False, index = False)

df = pd.read_csv(
    "all_sites_pieces/all_sites_03.vcf", sep = '\t', header = None, on_bad_lines = 'skip')
df = df.dropna(thresh = 143, axis = 0)
df.to_csv("all_sites_pieces/all_sites_03.vcf", sep = '\t', header = False, index = False)

df = pd.read_csv(
    "all_sites_pieces/all_sites_04.vcf", sep = '\t', header = None, on_bad_lines = 'skip')
df = df.dropna(thresh = 143, axis = 0)
df.to_csv("all_sites_pieces/all_sites_04.vcf", sep = '\t', header = False, index = False)

df = pd.read_csv(
    "all_sites_pieces/all_sites_05.vcf", sep = '\t', header = None, on_bad_lines = 'skip')
df = df.dropna(thresh = 143, axis = 0)
df.to_csv("all_sites_pieces/all_sites_05.vcf", sep = '\t', header = False, index = False)

df = pd.read_csv(
    "all_sites_pieces/all_sites_06.vcf", sep = '\t', header = None, on_bad_lines = 'skip')
df = df.dropna(thresh = 143, axis = 0)
df.to_csv("all_sites_pieces/all_sites_06.vcf", sep = '\t', header = False, index = False)

df = pd.read_csv(
    "all_sites_pieces/all_sites_07.vcf", sep = '\t', header = None, on_bad_lines = 'skip')
df = df.dropna(thresh = 143, axis = 0)
df.to_csv("all_sites_pieces/all_sites_07.vcf", sep = '\t', header = False, index = False)

df = pd.read_csv(
    "all_sites_pieces/all_sites_08.vcf", sep = '\t', header = None, on_bad_lines = 'skip')
df = df.dropna(thresh = 143, axis = 0)
df.to_csv("all_sites_pieces/all_sites_08.vcf", sep = '\t', header = False, index = False)

df = pd.read_csv(
    "all_sites_pieces/all_sites_09.vcf", sep = '\t', header = None, on_bad_lines = 'skip')
df = df.dropna(thresh = 143, axis = 0)
df.to_csv("all_sites_pieces/all_sites_09.vcf", sep = '\t', header = False, index = False)

df = pd.read_csv(
    "all_sites_pieces/all_sites_10.vcf", sep = '\t', header = None, on_bad_lines = 'skip')
df = df.dropna(thresh = 143, axis = 0)
df.to_csv("all_sites_pieces/all_sites_10.vcf", sep = '\t', header = False, index = False)

df = pd.read_csv(
    "all_sites_pieces/all_sites_11.vcf", sep = '\t', header = None, on_bad_lines = 'skip')
df = df.dropna(thresh = 143, axis = 0)
df.to_csv("all_sites_pieces/all_sites_11.vcf", sep = '\t', header = False, index = False)