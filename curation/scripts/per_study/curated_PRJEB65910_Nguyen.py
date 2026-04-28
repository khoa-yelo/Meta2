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


df_curated_PRJEB65910 = df_curated[df_curated['Study_ID'] == "PRJEB65910"].copy()

#create file
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB65910"][["sample_accession", "sample_sample-name"]]
df_original_subject["sample_sample-name"] = df_original_subject["sample_sample-name"].apply(lambda x: x)


sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_health_map = df_original[df_original['study_bioproject'] == "PRJEB65910"].set_index("sample_accession")["sample_sample-desc"].to_dict()
sample_accession_body_site_map = df_sra[df_sra['bioproject'] == "PRJNA801448"][["sample_accession", "body_site"]].set_index("sample_accession").to_dict()["body_site"]
sample_accession_location_map = df_sra[df_sra['bioproject'] == "PRJNA801448"][["sample_accession", "geo_location"]].set_index("sample_accession").to_dict()["geo_location"]
sample_accession_nucleotide_map = df_sra[df_sra['bioproject'] == "PRJNA801448"][["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]


df_curated_PRJEB65910["Subject_ID"] = df_curated_PRJEB65910["Sample_accession"].apply(lambda x: sample_accession_subject_map[x])
df_curated_PRJEB65910["Nucleotide_Type"] = df_curated_PRJEB65910["Sample_accession"].apply(lambda x: sample_accession_nucleotide_map[x])
df_curated_PRJEB65910["Body_site"] = df_curated_PRJEB65910["Sample_accession"].apply(lambda x: sample_accession_body_site_map[x])
df_curated_PRJEB65910["Body_site_core"] = "urinary"
df_curated_PRJEB65910["Health_status"] = df_curated_PRJEB65910["Sample_accession"].apply(lambda x: sample_accession_health_map[x])
df_curated_PRJEB65910["Location"] = df_curated_PRJEB65910["Sample_accession"].apply(lambda x: sample_accession_location_map[x])
df_curated_PRJEB65910.loc[df_curated_PRJEB65910["Subject_ID"] == "Water", "Age"] = np.nan
df_curated_PRJEB65910["Sex"] = "Female"
df_curated_PRJEB65910.loc[df_curated_PRJEB65910["Subject_ID"] == "Water", "Sex"] = np.nan
df_curated_PRJEB65910.loc[df_curated_PRJEB65910["Subject_ID"] == "Water", "Nucleotide_Type"] = np.nan

df_curated_PRJEB65910["Age_catagories"] = "Postmenopause"

df_curated_PRJEB65910.to_csv("./Final_Curated_Out/PRJEB65910_Nguyen.csv", index = None)
