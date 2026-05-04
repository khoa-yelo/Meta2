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

df_curated_PRJEB66073 = df_curated[df_curated['Study_ID'] == "PRJEB66073"]
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB66073"][["sample_accession", "sample_sample-name","sample_sample-desc","sample_biosample"]]
df_sra_subject = df_sra[df_sra['bioproject'] == "PRJNA752888"]
sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_nucleotide_map = df_sra_subject[["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]
sample_accession_age_map = ( df_sra_subject[["sample_accession", "age_years"]].set_index("sample_accession").to_dict()["age_years"])
sample_accession_sex_map = (df_sra_subject[["sample_accession", "sex"]].set_index("sample_accession").to_dict()["sex"])
sample_accession_diet_map = (df_sra_subject[["sample_accession", "diet"]].set_index("sample_accession").to_dict()["diet"])
sample_accession_health_map = (df_sra_subject[["sample_accession", "health_condition"]].set_index("sample_accession").to_dict()["health_condition"])
sample_accession_location_map = (df_sra_subject[["sample_accession", "geo_location"]].set_index("sample_accession").to_dict()["geo_location"])

df_curated_PRJEB66073["Subject_ID"] = df_curated_PRJEB66073["Sample_accession"].map(sample_accession_subject_map)
df_curated_PRJEB66073["Nucleotide_Type"] = df_curated_PRJEB66073["Sample_accession"].map(sample_accession_nucleotide_map)
df_curated_PRJEB66073["Body_site"] = "oral"
df_curated_PRJEB66073["Body_site_core"] = "oral"
#df_curated_PRJEB53689["Nucleotide_Type"] = "DNA"
df_curated_PRJEB66073["Age"] = df_curated_PRJEB66073["Sample_accession"].map(sample_accession_age_map)
df_curated_PRJEB66073["Sex"] = df_curated_PRJEB66073["Sample_accession"].map(sample_accession_sex_map)
df_curated_PRJEB66073["Health_status"] = df_curated_PRJEB66073["Sample_accession"].map(sample_accession_health_map)
df_curated_PRJEB66073["Location"] = df_curated_PRJEB66073['Sample_accession'].map(sample_accession_location_map)

df_curated_PRJEB66073.to_csv("/content/drive/Shareddrives/Meta2/Metadata/Final_Curated_Out/PRJEB66073.csv", index = None)
