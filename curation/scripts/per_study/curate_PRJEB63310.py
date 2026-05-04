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

df_curated_PRJEB63310 = df_curated[df_curated['Study_ID'] == "PRJEB63310"]
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB63310"][["sample_accession", "sample_sample-name","sample_sample-desc","sample_biosample"]]
df_sra_subject = df_sra[df_sra['bioproject'] == "PRJEB55713"]
sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_nucleotide_map = df_sra_subject[["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]
sample_accession_age_map = ( df_sra_subject[["sample_accession", "age_years"]].set_index("sample_accession").to_dict()["age_years"])
sample_accession_sex_map = (df_sra_subject[["sample_accession", "sex"]].set_index("sample_accession").to_dict()["sex"])
sample_accession_diet_map = (df_sra_subject[["sample_accession", "diet"]].set_index("sample_accession").to_dict()["diet"])
sample_accession_health_map = (df_sra_subject[["sample_accession", "health_condition"]].set_index("sample_accession").to_dict()["health_condition"])
sample_accession_location_map = (df_sra_subject[["sample_accession", "geo_location"]].set_index("sample_accession").to_dict()["geo_location"])
sample_accession_lifestyle_map = (df_sra_subject[["sample_accession", "lifestyle"]].set_index("sample_accession").to_dict()["lifestyle"])

df_sra[df_sra['bioproject'] == "PRJEB55713"] #use study_ID_original here

df_curated_PRJEB63310["Subject_ID"] = df_curated_PRJEB63310["Sample_accession"].map(sample_accession_subject_map)
df_curated_PRJEB63310["Nucleotide_Type"] = df_curated_PRJEB63310["Sample_accession"].map(sample_accession_nucleotide_map)
df_curated_PRJEB63310["Body_site"] = "gut"
df_curated_PRJEB63310["Body_site_core"] = "gut"
#df_curated_PRJEB53689["Nucleotide_Type"] = "DNA"
df_curated_PRJEB63310["Age"] = df_curated_PRJEB63310["Sample_accession"].map(sample_accession_age_map)
df_curated_PRJEB63310["Sex"] = df_curated_PRJEB63310["Sample_accession"].map(sample_accession_sex_map)
df_curated_PRJEB63310["Health_status"] = "Healthy"
df_curated_PRJEB63310["Location"] = df_curated_PRJEB63310['Sample_accession'].map(sample_accession_location_map)
df_curated_PRJEB63310.to_csv("/content/drive/Shareddrives/Meta2/Metadata/Final_Curated_Out/PRJEB63310.csv", index = None)
