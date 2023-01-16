
import pandas as pd
import re

# PXD007145 conversion
# no merge
data = pd.read_csv("D:/r-research/based-peptide-analysis/PXD007145/PXD007145-Th_filter.sdrf_openms_design_msstats_in"
                   ".csv",
                   sep=",", index_col=None)

data = data[-data["ProteinName"].str.contains("DECOY|CONTAMINANT")]


def sub_mod(peptide):
    peptide = peptide.replace(".", "")
    peptide = re.sub(r"\(.*?\)", "", peptide)
    return peptide


# data = data.groupby(["PeptideSequence", "ProteinName", "Reference", "Run"], as_index=False)["Intensity"].sum()

data.rename(columns={"ProteinName": "Proteins", "Intensity": "sumIntensity"}, inplace=True)

data["Run"] = data.apply(lambda x: "sumIntensity_" + str(x["Run"]), axis=1)

data = data.pivot(index=["Proteins", 'PeptideSequence', 'PrecursorCharge'], columns='Run', values='sumIntensity')
data = data.reset_index()
data["Sequence"] = data.apply(lambda x: sub_mod(x["PeptideSequence"]), axis=1)
data.to_csv("D:/r-research/PXD007145-Th_filter_no_merge.mzTab", sep="\t")

print(data.head())

# merge intensity of a peptide sequence with different charge

data = pd.read_csv("D:/r-research/based-peptide-analysis/PXD007145/PXD007145-Th_filter.sdrf_openms_design_msstats_in"
                   ".csv",
                   sep=",", index_col=None)

data = data[-data["ProteinName"].str.contains("DECOY|CONTAMINANT")]


def sub_mod(peptide):
    peptide = peptide.replace(".", "")
    peptide = re.sub(r"\(.*?\)", "", peptide)
    return peptide


data = data.groupby(["PeptideSequence", "ProteinName", "Reference", "Run"], as_index=False)["Intensity"].sum()

data.rename(columns={"ProteinName": "Proteins", "Intensity": "sumIntensity"}, inplace=True)

data["Run"] = data.apply(lambda x: "sumIntensity_" + str(x["Run"]), axis=1)

data = data.pivot(index=["Proteins", 'PeptideSequence'], columns='Run', values='sumIntensity')
data = data.reset_index()
data["Sequence"] = data.apply(lambda x: sub_mod(x["PeptideSequence"]), axis=1)
data.to_csv("D:/r-research/PXD007145-Th_filter_charge_merge.mzTab", sep="\t")

# merge intensity of peptide sequence with different charge and modification
data = pd.read_csv("D:/r-research/based-peptide-analysis/PXD007145/PXD007145-Th_filter.sdrf_openms_design_msstats_in"
                   ".csv",
                   sep=",", index_col=None)

data = data[-data["ProteinName"].str.contains("DECOY|CONTAMINANT")]


def sub_mod(peptide):
    peptide = peptide.replace(".", "")
    peptide = re.sub(r"\(.*?\)", "", peptide)
    return peptide


data["Sequence"] = data.apply(lambda x: sub_mod(x["PeptideSequence"]), axis=1)
data = data.groupby(["Sequence", "ProteinName", "Reference", "Run"], as_index=False)["Intensity"].sum()
data.rename(columns={"ProteinName": "Proteins", "Intensity": "sumIntensity"}, inplace=True)

data["Run"] = data.apply(lambda x: "sumIntensity_" + str(x["Run"]), axis=1)

data = data.pivot(index=["Proteins", 'Sequence'], columns='Run', values='sumIntensity')
data = data.reset_index()
data.to_csv("D:/r-research/PXD007145-Th_filter_all_merge.mzTab", sep="\t")


# PXD000279 conversion
data = pd.read_csv("C:/Users/ChengXin/Documents/GitHub/chengxin/quantms-research/r-research/based-peptide-analysis"
                   "/PXD000279/out_msstats_filter.csv",
                   sep=",", index_col=None)

data = data[-data["ProteinName"].str.contains("DECOY|CONTAMINANT")]


def sub_mod(peptide):
    peptide = peptide.replace(".", "")
    peptide = re.sub(r"\(.*?\)", "", peptide)
    return peptide


# data = data.groupby(["PeptideSequence", "ProteinName", "Reference", "Run"], as_index=False)["Intensity"].sum()

data.rename(columns={"ProteinName": "Proteins", "Intensity": "sumIntensity"}, inplace=True)

data["Run"] = data.apply(lambda x: "sumIntensity_" + str(x["Run"]), axis=1)

data = data.pivot(index=["Proteins", 'PeptideSequence', 'PrecursorCharge'], columns='Run', values='sumIntensity')
data = data.reset_index()
data["Sequence"] = data.apply(lambda x: sub_mod(x["PeptideSequence"]), axis=1)
data.to_csv("D:/r-research/PXD000279_filter_no_merge.sdrf_openms_design_msstats_in.csv", sep="\t")
