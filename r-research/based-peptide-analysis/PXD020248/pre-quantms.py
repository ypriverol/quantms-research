import re
import pandas as pd

out_msstats = pd.read_csv("./out_msstats.csv")
def sub_mod(peptide):
    peptide = peptide.replace(".", "")
    peptide = re.sub("\(.*?\)", "", peptide)
    return peptide
    
print(len(out_msstats))
out_msstats["sequence"] = out_msstats.apply(lambda x: sub_mod(x["PeptideSequence"]), axis=1)
filtered_df = out_msstats.groupby('sequence').filter(lambda x: len(set(x["ProteinName"])) == 1)
filtered_df = filtered_df[-filtered_df["ProteinName"].str.contains(";")]
filtered_df = filtered_df.groupby('ProteinName').filter(lambda x: len(set(x["sequence"])) > 2)
filtered_df.drop("sequence", axis=1, inplace=True)
print(len(filtered_df))
filtered_df.to_csv("./out_msstats_filter.csv", index=False)