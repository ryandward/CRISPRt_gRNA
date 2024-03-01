import pandas as pd
from find_crisprt_sites import find_crisprt_sites
from find_crisprt_sites import parse_transcription_units


df1 = find_crisprt_sites(
    "GCA_000005845.2.gb", "GCF_000005845.2_grna_top_ten_dedupe.tsv"
)

df2 = find_crisprt_sites(
    "GCA_003054575.1.gb", "GCF_003054575.1_grna_top_ten_dedupe.tsv"
)

df1_gene_to_tu = parse_transcription_units("EColi_TU.txt")
# somehow get the transcription units for ZMobilis

df1['TUTarget'] = df1['LocusTarget'].map(df1_gene_to_tu)
df1['TUInserted'] = df1['LocusInserted'].map(df1_gene_to_tu)
df1['SameTU'] = df1.apply(lambda row: None if pd.isnull(row['TUInserted']) else (True if row['TUTarget'] == row['TUInserted'] else False), axis=1)

# perform the same operations on df2
# df2['TUTarget'] = df2['LocusTarget'].map(df2_gene_to_tu)
# df2['TUInserted'] = df2['LocusInserted'].map(df2_gene_to_tu)
# df2['SameTU'] = df2.apply(lambda row: None if pd.isnull(row['TUInserted']) else (True if row['TUTarget'] == row['TUInserted'] else False), axis=1)

df1.to_csv("EColi_CtSites.tsv", index=False, sep="\t")
df2.to_csv("ZMobilis_CtSites.tsv", index=False, sep="\t")
