import numpy as np
import pandas as pd
from GenBankParser import GenBankParser
from find_crisprt_sites import find_crisprt_sites
from find_crisprt_sites import parse_transcription_units


df1 = find_crisprt_sites(
    "GCA_000005845.2.gb", "GCF_000005845.2_grna_top_ten_dedupe.tsv"
)

df1_gene_to_tu = parse_transcription_units("EColi_TU.txt")
# somehow get the transcription units for ZMobilis

df1["TUTarget"] = df1["LocusTarget"].map(df1_gene_to_tu)
df1["TUInserted"] = df1["LocusInserted"].map(df1_gene_to_tu)
df1["SameTU"] = df1.apply(
    lambda row: (
        None
        if pd.isnull(row["TUInserted"])
        else (True if row["TUTarget"] == row["TUInserted"] else False)
    ),
    axis=1,
)

# extract the start site for each "LocusTarget", then find the distance between "InsertSite" and "LocusTargetStart
gbp1 = GenBankParser("GCA_000005845.2.gb")
# Prepare mappings for both "Start" and "End" indexed by "Locus_Tag"
df1_start_mapping = gbp1.ranges.df.loc[gbp1.ranges.df["Type"] == "gene"].set_index(
    "Locus_Tag"
)["Start"]
df1_end_mapping = gbp1.ranges.df.loc[gbp1.ranges.df["Type"] == "gene"].set_index(
    "Locus_Tag"
)["End"]

# Apply mapping with validation
df1["LocusTargetStart"] = np.where(
    df1["CtStrand"] == "+",
    df1["LocusTarget"].map(df1_start_mapping),
    np.where(df1["CtStrand"] == "-", df1["LocusTarget"].map(df1_end_mapping), 0),
)

# Calculate "Distance" with validation
df1["CtUpstream"] = np.where(
    df1["CtStrand"] == "+",
    df1["LocusTargetStart"] - df1["InsertSite"],
    np.where(
        df1["CtStrand"] == "-",
        df1["InsertSite"] - df1["LocusTargetStart"],
        0,
    ),
)

df1.to_csv("EColi_CtSites.tsv", index=False, sep="\t")

df2 = find_crisprt_sites(
    "GCA_003054575.1.gb", "GCF_003054575.1_grna_top_ten_dedupe.tsv"
)
# perform the same operations on df2
# df2['TUTarget'] = df2['LocusTarget'].map(df2_gene_to_tu)
# df2['TUInserted'] = df2['LocusInserted'].map(df2_gene_to_tu)
# df2['SameTU'] = df2.apply(lambda row: None if pd.isnull(row['TUInserted']) else (True if row['TUTarget'] == row['TUInserted'] else False), axis=1)
gbp2 = GenBankParser("GCA_003054575.1.gb")

zmo_last_compound_location_subset = gbp2.ranges.df.loc[
    gbp2.ranges.df["Type"] == "gene"
].drop_duplicates(subset="Locus_Tag", keep="last")


# Prepare mappings for both "Start" and "End" indexed by "Locus_Tag"
df2_start_mapping = zmo_last_compound_location_subset.loc[
    zmo_last_compound_location_subset["Type"] == "gene"
].set_index("Locus_Tag")["Start"]
df2_end_mapping = zmo_last_compound_location_subset.loc[
    zmo_last_compound_location_subset["Type"] == "gene"
].set_index("Locus_Tag")["End"]

# Apply mapping with validation
df2["LocusTargetStart"] = np.where(
    df2["CtStrand"] == "+",
    df2["LocusTarget"].map(df2_start_mapping),
    np.where(df2["CtStrand"] == "-", df2["LocusTarget"].map(df2_end_mapping), 0),
)

# Calculate "Distance" with validation
df2["CtUpstream"] = np.where(
    df2["CtStrand"] == "+",
    df2["LocusTargetStart"] - df2["InsertSite"],
    np.where(
        df2["CtStrand"] == "-",
        df2["InsertSite"] - df2["LocusTargetStart"],
        0,
    ),
)

df2.to_csv("ZMobilis_CtSites.tsv", index=False, sep="\t")
