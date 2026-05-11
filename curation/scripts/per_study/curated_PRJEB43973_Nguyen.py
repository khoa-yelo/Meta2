# mount drive
from google.colab import drive
drive.mount('/content/drive')


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


sra_cleaned_path = "/content/drive/Shareddrives/Meta2/Metadata/Curated_org_data/cleaned_sra_metadata.tsv"
curated_cleaned_path = "/content/drive/Shareddrives/Meta2/Metadata/Curated_org_data/curated_metadata.csv"
original_metadata_path = "/content/drive/Shareddrives/Meta2/Data/hmb_assemblies_metadata.csv"

df_sra = pd.read_csv(sra_cleaned_path, sep="\t", low_memory=False)
df_curated = pd.read_csv(curated_cleaned_path, sep=",")
df_original = pd.read_csv(original_metadata_path, sep=",", low_memory=False)


df_curated_PRJEB43973 = df_curated[df_curated['Study_ID'] == "PRJEB43973"].copy()

#create file
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB43973"][["sample_accession", "sample_sample-name"]]
df_original_subject["sample_sample-name"] = df_original_subject["sample_sample-name"].apply(lambda x: x)

sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_body_site_map = df_sra[df_sra['bioproject'] == "PRJNA385949"][["sample_accession", "body_site"]].set_index("sample_accession").to_dict()["body_site"]
sample_accession_location_map = df_sra[df_sra['bioproject'] == "PRJNA385949"][["sample_accession", "geo_location"]].set_index("sample_accession").to_dict()["geo_location"]
sample_accession_nucleotide_map = df_sra[df_sra['bioproject'] == "PRJNA385949"][["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]
sample_accession_health_map = df_original[df_original['study_bioproject'] == "PRJEB43973"].set_index("sample_accession")["sample_sample-desc"].to_dict()

df_curated_PRJEB43973["Subject_ID"] = df_curated_PRJEB43973["Sample_accession"].apply(lambda x: sample_accession_subject_map[x])
df_curated_PRJEB43973["Nucleotide_Type"] = df_curated_PRJEB43973["Sample_accession"].apply(lambda x: sample_accession_nucleotide_map[x])
df_curated_PRJEB43973["Body_site"] = df_curated_PRJEB43973["Sample_accession"].apply(lambda x: sample_accession_body_site_map[x])
df_curated_PRJEB43973["Body_site_core"] = "Gut"
df_curated_PRJEB43973["Location"] = df_curated_PRJEB43973["Sample_accession"].apply(lambda x: sample_accession_location_map[x])
df_curated_PRJEB43973["Health_status"] = df_curated_PRJEB43973["Sample_accession"].map(sample_accession_health_map)


df_curated_PRJEB43973["Age_catagories"] = np.nan
df_curated_PRJEB43973["Age_catagories"] = df_curated_PRJEB43973["Age_catagories"].astype(object)

df_curated_PRJEB43973.loc[df_curated_PRJEB43973["Subject_ID"].str.startswith("p", na=False), "Age_catagories"] = "from young_adults to middle_age"
df_curated_PRJEB43973.loc[df_curated_PRJEB43973["Subject_ID"].str.startswith("SKST", na=False), "Age_catagories"] = "children"
df_curated_PRJEB43973["Age_catagories"] = df_curated_PRJEB43973["Age_catagories"].astype(object)

df_curated_PRJEB43973.to_csv("./Final_Curated_Out/PRJEB43973_Nguyen.csv", index = None)

