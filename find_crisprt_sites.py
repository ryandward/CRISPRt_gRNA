from GenBankParser import GenBankParser
import pandas as pd
import pyranges as pr


def find_crisprt_sites(genbank_file: str, crisprt_file: str) -> pd.DataFrame:
    gbp = GenBankParser(genbank_file)
    CtTargets = pd.read_csv(crisprt_file, sep="\t")

    CtTargets = CtTargets.rename(
        columns={
            "chrom": "Chromosome",
            "start": "Start",
            "end": "End",
            "gene": "LocusTarget",
            "target": "CtTarget",
        }
    )

    CtTargets["GeneTarget"] = CtTargets["LocusTarget"].apply(
        gbp.find_gene_name_for_locus
    )
    

    CtTargets["Strand"] = CtTargets["repldir"].apply(
        lambda x: "+" if x == "fwd" else "-"
    )

    CtTargets = CtTargets.astype({"Start": int, "End": int})

    CtTargets["CtStart"] = CtTargets["Start"].astype(str)
    CtTargets["CtEnd"] = CtTargets["End"].astype(str)

    for index, row in CtTargets.iterrows():
        if row["Strand"] == "+":
            CtTargets.at[index, "Start"] = row["End"] + 49
            CtTargets.at[index, "End"] = CtTargets.at[index, "Start"]
        elif row["Strand"] == "-":
            CtTargets.at[index, "Start"] = row["Start"] - 49
            CtTargets.at[index, "End"] = CtTargets.at[index, "Start"]
        else:
            raise ValueError("Strand must be either + or -")

    CtTargets = CtTargets[
        [
            "CtTarget",
            "LocusTarget",
            "GeneTarget",
            "Chromosome",
            "Start",
            "End",
            "Strand",
            "CtStart",
            "CtEnd",
        ]
    ]
    CtTargets = CtTargets.drop_duplicates()
    CtTargets = CtTargets.dropna()

    CtTargets = pr.PyRanges(CtTargets)

    gbp.geneRanges = gbp.ranges[gbp.ranges.Type == "gene"]

    sites = CtTargets.join(gbp.geneRanges, how="left", suffix="Inserted")
    sites = sites.sort()
    sites.InsertSite = sites.Start
    sites = sites[
        [
            "CtTarget",
            "LocusTarget",
            "GeneTarget",
            "Chromosome",
            "CtStart",
            "CtEnd",
            "Strand",
            "InsertSite",
            "Locus_Tag",
            "Gene",
            "StrandInserted",
        ]
    ]

    df = sites.df
    df = df.rename(
        columns={
            "Locus_Tag": "LocusInserted",
            "Gene": "GeneInserted",
            "StrandInserted": "StrandInserted",
            "Strand": "CtStrand",
        }
    )

    # Drop Start and End columns
    df = df.drop(columns=["Start", "End"])

    # if LocusInserted is -1, make it NaN
    df["LocusInserted"] = df["LocusInserted"].apply(lambda x: None if x == "-1" else x)

    df["GeneInserted"] = df["GeneInserted"].apply(lambda x: None if x == "-1" else x)

    return df

def parse_transcription_units(file_path):
    gene_to_tu = {}
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            tu_id = fields[0]
            if len(fields) > 1:
                genes = fields[1].split(' // ')
                for gene in genes:
                    gene_to_tu[gene] = tu_id
    return gene_to_tu

