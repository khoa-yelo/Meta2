import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


sra_cleaned_path = "./Curated_org_data/cleaned_sra_metadata.tsv"
curated_cleaned_path = "./Curated_org_data/curated_metadata.csv"
original_metadata_path = "./Curated_org_data/hmb_assemblies_metadata.csv"

df_sra = pd.read_csv(sra_cleaned_path, sep="\t", low_memory=False)
df_curated = pd.read_csv(curated_cleaned_path, sep=",")
df_original = pd.read_csv(original_metadata_path, sep=",", low_memory=False)


df_curated_PRJEB65450 = df_curated[df_curated['Study_ID'] == "PRJEB65450"].copy()

#create file
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB65450"][["sample_accession", "sample_sample-name"]]
df_original_subject["sample_sample-name"] = df_original_subject["sample_sample-name"].apply(lambda x: x)

sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_sex_map = df_original[df_original['study_bioproject'] == "PRJEB65450"].set_index("sample_accession")["sample_sample-desc"].str.extract(r'(?i)(male|female)', expand=False).str.capitalize().to_dict()
sample_accession_location_map = df_sra[df_sra['bioproject'] == "PRJEB54966"][["sample_accession", "geo_location"]].set_index("sample_accession").to_dict()["geo_location"]
sample_accession_nucleotide_map = df_sra[df_sra['bioproject'] == "PRJEB54966"][["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]

df_curated_PRJEB65450["Subject_ID"] = df_curated_PRJEB65450["Sample_accession"].apply(lambda x: sample_accession_subject_map[x])
df_curated_PRJEB65450["Nucleotide_Type"] = df_curated_PRJEB65450["Sample_accession"].apply(lambda x: sample_accession_nucleotide_map[x])
df_curated_PRJEB65450["Location"] = df_curated_PRJEB65450["Sample_accession"].apply(lambda x: sample_accession_location_map[x])
df_curated_PRJEB65450["Sex"] = df_curated_PRJEB65450["Sample_accession"].apply(lambda x: sample_accession_sex_map[x])

df_curated_PRJEB65450.to_csv("./Final_Curated_Out/PRJEB65450_Nguyen.csv", index = None)
