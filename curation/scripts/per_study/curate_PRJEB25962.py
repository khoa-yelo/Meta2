 # mount drive
from google.colab import drive
drive.mount('/content/drive')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sra_cleaned_path = "/content/drive/Shareddrives/Meta2/Metadata/cleaned_sra_metadata.tsv"
curated_cleaned_path = "/content/drive/Shareddrives/Meta2/Metadata/curated_metadata.csv"
original_metadata_path = "/content/drive/Shareddrives/Meta2/Data/hmb_assemblies_metadata.csv"

df_sra = pd.read_csv(sra_cleaned_path, sep="\t")
df_curated = pd.read_csv(curated_cleaned_path, sep=",")
df_original = pd.read_csv(original_metadata_path, sep=",")

df_curated_PRJEB25962 = df_curated[df_curated['Study_ID'] == "PRJEB25962"]
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB25962"][["sample_accession", "sample_sample-name","sample_sample-desc","sample_biosample"]]
df_sra_subject = df_sra[df_sra['bioproject'] == "PRJEB1220"]
sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_nucleotide_map = df_sra_subject[["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]
sample_accession_age_map = ( df_sra_subject[["sample_accession", "Age"]].set_index("sample_accession").to_dict()["Age"])
sample_accession_sex_map = (df_sra_subject[["sample_accession", "Gender"]].set_index("sample_accession").to_dict()["Gender"])
sample_accession_health_map = (df_sra_subject[["sample_accession", "Health Status"]].set_index("sample_accession").to_dict()["Health Status"])
sample_accession_location_map = (df_sra_subject[["sample_accession", "Nationality"]].set_index("sample_accession").to_dict()["Nationality"])

df_sra[df_sra['bioproject'] == "PRJEB1220"] #use study_ID_original here

df_original_subject['sample_sample-name'].unique()

PRJEB25962_meta = pd.read_excel("/content/drive/Shareddrives/Meta2/Metadata/Manual_metadata/PRJEB25962.xls", header = None)

df_original_subject['sample_sample-name'] = df_original_subject['sample_sample-name'].str.replace('-', '.', regex=False)
df_original_subject

new_header = PRJEB25962_meta.iloc[0] 
PRJEB25962_meta = PRJEB25962_meta[1:] 
PRJEB25962_meta.columns = new_header
PRJEB25962_meta.reset_index(drop=True, inplace=True)

PRJEB25962_meta = PRJEB25962_meta.rename(columns={"Sample ID" : 'sample_sample-name', "Individual ID" : "Individual_ID"})
df_with_meta = df_original_subject.merge(PRJEB25962_meta, on="sample_sample-name", how="left")

df_with_meta['Gender'].replace({'F': 'Female', 'M':'Male'}, inplace=True)
df_with_meta['Nationality'] = df_with_meta['Nationality'].str.capitalize()
mapping = {
    'Spanish': 'Spain',
    'Danish': 'Denmark'
}

df_with_meta['Nationality'] = df_with_meta['Nationality'].replace(mapping)
df_with_meta


df_sra_subject = df_sra_subject.merge(
    df_with_meta[["sample_accession", "Age", "Gender", "Nationality", "Health Status"]],
    on="sample_accession",
    how="left"
)
df_curated_PRJEB25962["Subject_ID"] = df_curated_PRJEB25962["Sample_accession"].map(sample_accession_subject_map)
df_curated_PRJEB25962["Nucleotide_Type"] = df_curated_PRJEB25962["Sample_accession"].map(sample_accession_nucleotide_map)
df_curated_PRJEB25962["Body_site"] = "gut"
df_curated_PRJEB25962["Body_site_core"] = "gut"
#df_curated_PRJEB53689["Nucleotide_Type"] = "DNA"
df_curated_PRJEB25962["Age"] = df_curated_PRJEB25962["Sample_accession"].map(sample_accession_age_map)
df_curated_PRJEB25962["Sex"] = df_curated_PRJEB25962["Sample_accession"].map(sample_accession_sex_map)
df_curated_PRJEB25962["Health_status"] = df_curated_PRJEB25962["Sample_accession"].map(sample_accession_health_map)
df_curated_PRJEB25962["Location"] = df_curated_PRJEB25962['Sample_accession'].map(sample_accession_location_map)

df_curated_PRJEB25962.to_csv("/content/drive/Shareddrives/Meta2/Metadata/Final_Curated_Out/PRJEB25962.csv", index = None)
