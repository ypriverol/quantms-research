import re
import pandas as pd
from collections import Counter

def clean_PeptideSequence(pep_seq):
    pattern = re.compile(r"\(.*?\)")
    pep_seq = pep_seq.replace(".", "")
    return pattern.sub("", pep_seq)

def modified_PeptideSequence(pep_seq):
    pattern = re.compile(r"\(.*?\)")
    t_pep_seq = pep_seq.replace(".", "").replace("Acetyl", "ac").replace("Carbamidomethyl", "ca").replace("Oxidation", "ox")
    t_pep_seq = "_" + t_pep_seq + "_"
    return t_pep_seq

def get_modifications(pep_seq):
    pattern = re.compile(r"\(.*?\)")
    t = ""
    modifications = re.findall(pattern, pep_seq)
    if len(modifications) == 0:
        return "Unmodified"
    c = Counter(modifications)
    for key, value in c.items():
        if value == 1:
            c[key] = ""
        else:
            c[key] = str(c[key]) + " "
    modifications = ",".join(modifications)
    if "(Acetyl)" in c:
        modifications = modifications.replace("(Acetyl)", c["(Acetyl)"] + "Acetyl (Protein N-term)")
    if "(Carbamidomethyl)" in c:
        modifications = modifications.replace("(Carbamidomethyl)", c["(Carbamidomethyl)"] + "Carbamidomethyl (C)")
    if "(Oxidation)" in c:
        modifications = modifications.replace("(Oxidation)", c["(Oxidation)"] + "Oxidation (M)")
        
    return modifications

def get_experiment(row):
    return "_".join(row['Reference'].split("_")[-3:-1])
    
    
evi = pd.DataFrame(None, columns=['sequence', 'modified_sequence', 'modifications', 'protein_group', 
                                  'protein', 'experiment', 'charge', 'reverse', 'contaminant', 'intensity'])
quantms = pd.read_csv("./out_msstats_filter.csv", sep=',', header=0)
quantms = quantms[-quantms['ProteinName'].str.contains("DECOY_")]
evi['sequence'] = quantms.apply(lambda x: clean_PeptideSequence(x['PeptideSequence']), axis=1)
evi['modified_sequence'] = quantms.apply(lambda x: modified_PeptideSequence(x['PeptideSequence']), axis=1)
evi['modifications'] = quantms.apply(lambda x: get_modifications(x['PeptideSequence']), axis=1)
evi['protein_group'] = quantms['ProteinName']
evi['protein'] = quantms['ProteinName']
evi['experiment'] = quantms.apply(lambda x: get_experiment(x), axis=1)
evi['charge'] = quantms['PrecursorCharge']
evi['intensity'] = quantms['Intensity']
evi.to_csv("./out_proteus.csv", sep=',', index=False)