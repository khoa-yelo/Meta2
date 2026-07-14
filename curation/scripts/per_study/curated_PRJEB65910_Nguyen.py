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

df_sra = pd.read_csv(sra_cleaned_path, sep="\t")
df_curated = pd.read_csv(curated_cleaned_path, sep=",")
df_original = pd.read_csv(original_metadata_path, sep=",")


df_curated_PRJEB65910 = df_curated[df_curated['Study_ID'] == "PRJEB65910"].copy()

#check
print(df_curated.columns)
for col in df_curated[df_curated["Study_ID"] == "PRJEB65910"].columns:
  print(df_curated[df_curated["Study_ID"] == "PRJEB65910"][col].value_counts())
print(df_curated[df_curated["Study_ID"] == "PRJEB65910"]["Body_site"].value_counts(dropna=False))
print(df_curated[df_curated["Study_ID"] == "PRJEB65910"]["Body_site_core"].value_counts(dropna=False))
print(df_curated[df_curated["Study_ID"] == "PRJEB65910"]["Health_status"].value_counts(dropna=False))
print(df_curated[df_curated["Study_ID"] == "PRJEB65910"]["Diet"].value_counts(dropna=False))
print(df_curated[df_curated["Study_ID"] == "PRJEB65910"]["Location"].value_counts(dropna=False))
print(df_curated[df_curated["Study_ID"] == "PRJEB65910"]["Lifestyle"].value_counts(dropna=False))
print(df_curated[df_curated["Study_ID"] == "PRJEB65910"]["Sex"].value_counts(dropna=False))
print(df_curated[df_curated["Study_ID"] == "PRJEB65910"]["Age"].value_counts(dropna=False))
print(df_curated[df_curated["Study_ID"] == "PRJEB65910"]["Nucleotide_Type"].value_counts(dropna=False))

print(df_original.columns)
for col in df_original.columns:
    print(col)
pd.set_option('display.max_rows', 500)
subset = df_original[df_original["study_bioproject"] == "PRJEB65910"]
for col in subset.columns:
    counts = subset[col].value_counts(dropna=True)
    if not counts.empty:
        print(f"\n{'='*10} {col} {'='*10}")
        if len(counts) > 20:
            print(counts.head(20))
            print(f"... Total {len(counts)} unique values")
        else:
            print(counts)
df_original[df_original["study_bioproject"] == "PRJEB65910"]["sample_samplemeta_host_sex"]
print(df_original[df_original["study_bioproject"] == "PRJEB65910"]["sample_samplemeta_host_sex"].value_counts(dropna=False))
print(df_original[df_original["study_bioproject"] == "PRJEB65910"]["sample_samplemeta_body_site"].value_counts(dropna=False))
print(df_original[df_original["study_bioproject"] == "PRJEB65910"]["sample_samplemeta_age"].value_counts(dropna=False))
print(df_original[df_original["study_bioproject"] == "PRJEB65910"]["sample_samplemeta_diet"].value_counts(dropna=False))
print(df_original[df_original["study_bioproject"] == "PRJEB65910"]["sample_samplemeta_nucleic_acid_extraction"].value_counts(dropna=False))

df_sra.columns
for col in df_sra[df_sra['bioproject'] == "PRJNA801448"].columns:
  print(df_sra[df_sra['bioproject'] == "PRJNA801448"][col].value_counts())
df_sra[df_sra['bioproject'] == "PRJNA801448"]["biosample"].value_counts(dropna=False)
print(df_sra[df_sra['bioproject'] == "PRJNA801448"]["body_site"].value_counts(dropna=False))
df_sra[df_sra['bioproject'] == "PRJNA801448"]["age_years"].value_counts(dropna=False)
df_sra[df_sra['bioproject'] == "PRJNA801448"]["sex"].value_counts(dropna=False)
df_sra[df_sra['bioproject'] == "PRJNA801448"]["health_condition_potential"].value_counts(dropna=False)
df_sra[df_sra['bioproject'] == "PRJNA801448"]["diet"].value_counts(dropna=False)
df_sra[df_sra['bioproject'] == "PRJNA801448"]["geo_location"].value_counts(dropna=False)
df_sra[df_sra['bioproject'] == "PRJNA801448"]["lifestyle"].value_counts(dropna=False)
df_sra[df_sra['bioproject'] == "PRJNA801448"]["nucleotide_type"].value_counts(dropna=False)

#create file
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB65910"][["sample_accession", "sample_sample-name"]]
df_original_subject["sample_sample-name"] = df_original_subject["sample_sample-name"].apply(lambda x: x)

sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_body_site_map = df_sra[df_sra['bioproject'] == "PRJNA801448"][["sample_accession", "body_site"]].set_index("sample_accession").to_dict()["body_site"]
sample_accession_health_condition_potential_map = df_sra[df_sra['bioproject'] == "PRJNA801448"][["sample_accession", "health_condition_potential"]].set_index("sample_accession").to_dict()["health_condition_potential"]
sample_accession_location_map = df_sra[df_sra['bioproject'] == "PRJNA801448"][["sample_accession", "geo_location"]].set_index("sample_accession").to_dict()["geo_location"]
sample_accession_nucleotide_map = df_sra[df_sra['bioproject'] == "PRJNA801448"][["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]


df_curated_PRJEB65910["Subject_ID"] = df_curated_PRJEB65910["Sample_accession"].apply(lambda x: sample_accession_subject_map[x])
df_curated_PRJEB65910["Nucleotide_Type"] = df_curated_PRJEB65910["Sample_accession"].apply(lambda x: sample_accession_nucleotide_map[x])
df_curated_PRJEB65910["Body_site"] = df_curated_PRJEB65910["Sample_accession"].apply(lambda x: sample_accession_body_site_map[x])
df_curated_PRJEB65910["Health_status"] = df_curated_PRJEB65910["Sample_accession"].apply(lambda x: sample_accession_health_condition_potential_map[x])
df_curated_PRJEB65910["Location"] = df_curated_PRJEB65910["Sample_accession"].apply(lambda x: sample_accession_location_map[x])
df_curated_PRJEB65910["Age"] = "Postmenopause"
df_curated_PRJEB65910.loc[df_curated_PRJEB65910["Subject_ID"] == "Water", "Age"] = np.nan
df_curated_PRJEB65910["Sex"] = "Female"
df_curated_PRJEB65910.loc[df_curated_PRJEB65910["Subject_ID"] == "Water", "Sex"] = np.nan
df_curated_PRJEB65910.loc[df_curated_PRJEB65910["Subject_ID"] == "Water", "Nucleotide_Type"] = np.nan

df_curated_PRJEB65910.to_csv("/content/drive/Shareddrives/Meta2/Metadata/Final_Curated_Out/PRJEB65910_Nguyen.csv", index = None)
