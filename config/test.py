import pandas as pd
df = pd.read_csv("samples.tsv", sep="\t")
print(df.columns)
print(df.query("ragtag == 'Y'"))
