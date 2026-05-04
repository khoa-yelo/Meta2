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

df_curated_PRJEB65452 =  df_curated[df_curated['Study_ID'] == "PRJEB65452"]
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB65452"][["sample_accession", "sample_sample-name","sample_sample-desc","sample_biosample"]]
df_sra_subject = df_sra[df_sra['bioproject'] == "PRJNA932553"]
sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_nucleotide_map = df_sra_subject[["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]
sample_accession_age_map = ( df_sra_subject[["sample_accession", "age_years"]].set_index("sample_accession").to_dict()["Age"])
sample_accession_sex_map = (df_sra_subject[["sample_accession", "sex"]].set_index("sample_accession").to_dict()["Sex"])
sample_accession_health_map = (df_sra_subject[["sample_accession", "Status"]].set_index("sample_accession").to_dict()["Status"])

df_sra[df_sra['bioproject'] == "PRJNA932553"] #use study_ID_original here
df_original_subject['Patient_ID'] = ['C_3', 'P_9', 'P_2', 'P_2', 'P_15', 'P_14', 'P_20', 'C_1', 'C_5','P_16', 'P_20','C_5', 'P_5', 'P_1', 'P_19',
                                     'P_5', 'P_23', 'C_9', 'P_10', 'C_10', 'P_21', 'P_22', 'C_10', 'C_4', 'P_18', 'P_7', 'P_12', 'C_15', 'C_13',
                                     'C_6', 'P_3', 'C_11', 'C_14', 'P_11', 'P_8', 'P_23', 'C_12', 'C_2', 'P_8', 'C_16', 'P_9', 'P_4', 'P_13', 'C_7',
                                     'P_17', 'P_4', 'P_6', 'C_8']
PRJEB65452_meta = pd.read_excel("/content/drive/Shareddrives/Meta2/Metadata/Manual_metadata/PRJEB65452.xlsx", header = None)
PRJEB65452_meta = PRJEB65452_meta.iloc[2:]
new_header = PRJEB65452_meta.iloc[0] 
PRJEB65452_meta = PRJEB65452_meta[1:] 
PRJEB65452_meta.columns = new_header 
PRJEB65452_meta.reset_index(drop=True, inplace=True)
PRJEB65452_meta = PRJEB65452_meta.rename(columns={'Patient ID': 'Patient_ID'})
PRJEB65452_meta['Sex'] = PRJEB65452_meta['Sex'].map({'F': 'Female', 'M':'Male'})
PRJEB65452_meta['Status'] = PRJEB65452_meta['Status'].map({'Perio': 'Periodontitis', 'Control':'Control'})

df_with_meta = df_original_subject.merge(PRJEB65452_meta, on="Patient_ID", how="left")
df_sra_subject = df_sra_subject.merge(
    df_with_meta[["sample_accession", "Age", "Sex", "Status"]],
    on="sample_accession",
    how="left"
)

df_curated_PRJEB65452["Subject_ID"] = df_curated_PRJEB65452["Sample_accession"].map(sample_accession_subject_map)
df_curated_PRJEB65452["Nucleotide_Type"] = df_curated_PRJEB65452["Sample_accession"].map(sample_accession_nucleotide_map)
df_curated_PRJEB65452["Location"] = "Hong Kong"
df_curated_PRJEB65452["Body_site"] = "oral"
df_curated_PRJEB65452["Body_site_core"] = "oral"
df_curated_PRJEB65452["Age"] = df_curated_PRJEB65452["Sample_accession"].map(sample_accession_age_map)
df_curated_PRJEB65452["Sex"] = df_curated_PRJEB65452["Sample_accession"].map(sample_accession_sex_map)
df_curated_PRJEB65452["Health_status"] = df_curated_PRJEB65452["Sample_accession"].map(sample_accession_health_map)

df_curated_PRJEB65452.to_csv("/content/drive/Shareddrives/Meta2/Metadata/Final_Curated_Out/PRJEB65452.csv", index = None)
