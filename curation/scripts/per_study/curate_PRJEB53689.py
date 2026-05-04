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

df_curated_PRJEB53689 =  df_curated[df_curated['Study_ID'] == "PRJEB53689"]
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB53689"][["sample_accession", "sample_sample-name","sample_sample-desc","sample_biosample"]]
df_sra_subject = df_sra[df_sra['bioproject'] == "PRJEB51898"]
sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_nucleotide_map = df_sra_subject[["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]
sample_accession_age_map = ( df_sra_subject[["sample_accession", "Age"]].set_index("sample_accession").to_dict()["Age"])
sample_accession_sex_map = (df_sra_subject[["sample_accession", "Gender"]].set_index("sample_accession").to_dict()["Gender"])
sample_accession_diet_map = (df_sra_subject[["sample_accession", "Diet"]].set_index("sample_accession").to_dict()["Diet"])

PRJEB53689_meta = pd.read_excel("/content/drive/Shareddrives/Meta2/Metadata/Manual_metadata/PRJEB53689.xlsx", header = None)
PRJEB53689_meta = PRJEB53689_meta.iloc[2:]
new_header = PRJEB53689_meta.iloc[0] 
PRJEB53689_meta = PRJEB53689_meta[1:] 
PRJEB53689_meta.columns = new_header 
PRJEB53689_meta.reset_index(drop=True, inplace=True)
PRJEB53689_meta = PRJEB53689_meta.rename(columns={"Subject ID" : 'Subject_ID'})
PRJEB53689_meta = PRJEB53689_meta.rename(columns={"Dietary mode" : 'Diet'})

df_original_subject['Subject_ID'] = (
    df_original_subject['sample_sample-name']
    .str.replace(r'([GP])-vig-', r'\1V', regex=True)
)
df_original_subject
df_with_meta = df_original_subject.merge(PRJEB53689_meta, on="Subject_ID", how="left")


df_sra_subject = df_sra_subject.merge(
    df_with_meta[["sample_accession", "Age", "Gender", "Diet"]],
    on="sample_accession",
    how="left"
)





df_curated_PRJEB53689["Subject_ID"] = df_curated_PRJEB53689["Sample_accession"].map(sample_accession_subject_map)
df_curated_PRJEB53689["Nucleotide_Type"] = df_curated_PRJEB53689["Sample_accession"].map(sample_accession_nucleotide_map)
df_curated_PRJEB53689["Location"] = "China"
df_curated_PRJEB53689["Body_site"] = "vaginal"
df_curated_PRJEB53689["Body_site_core"] = "vaginal"
df_curated_PRJEB53689["Nucleotide_Type"] = "DNA"
df_curated_PRJEB53689["Health_status"] = "Healthy"
df_curated_PRJEB53689["Lifestyle"] = "no medical conditions, no antibiotic or microbial modulator within 3 months before sampling, abstinance and no vaginal douche for 5 days prior to sampling"


df_curated_PRJEB53689["Subject_ID"] = df_curated_PRJEB53689["Sample_accession"].map(sample_accession_subject_map)
df_curated_PRJEB53689["Nucleotide_Type"] = df_curated_PRJEB53689["Sample_accession"].map(sample_accession_nucleotide_map)
df_curated_PRJEB53689["Body_site"] = "vaginal"
df_curated_PRJEB53689["Body_site_core"] = "vaginal"
df_curated_PRJEB53689["Age"] = df_curated_PRJEB53689["Sample_accession"].map(sample_accession_age_map)
df_curated_PRJEB53689["Sex"] = df_curated_PRJEB53689["Sample_accession"].map(sample_accession_sex_map)
df_curated_PRJEB53689["Diet"] = df_curated_PRJEB53689["Sample_accession"].map(sample_accession_diet_map)
df_curated_PRJEB53689["Location"] = "China"

df_curated_PRJEB53689.to_csv("/content/drive/Shareddrives/Meta2/Metadata/Final_Curated_Out/PRJEB53689.csv", index = None)
