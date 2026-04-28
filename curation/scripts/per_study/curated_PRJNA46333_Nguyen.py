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


df_curated_PRJEB33701 = df_curated[df_curated['Study_ID'] == "PRJEB33701"].copy()

#create file
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB33701"][["sample_accession", "sample_sample-name"]]
df_original_subject["sample_sample-name"] = df_original_subject["sample_sample-name"].apply(lambda x: x)

sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_body_site_map = df_sra[df_sra['bioproject'] == "PRJNA46333"][["sample_accession", "body_site"]].set_index("sample_accession").to_dict()["body_site"]
sample_accession_sex_map = df_sra[df_sra['bioproject'] == "PRJNA46333"][["sample_accession", "sex"]].set_index("sample_accession").to_dict()["sex"]
sample_accession_location_map = df_sra[df_sra['bioproject'] == "PRJNA46333"][["sample_accession", "geo_location"]].set_index("sample_accession").to_dict()["geo_location"]
sample_accession_health_map = df_sra[df_sra['bioproject'] == "PRJNA46333"][["sample_accession", "health_condition_potential"]].set_index("sample_accession").to_dict()["health_condition_potential"]
sample_accession_nucleotide_map = df_sra[df_sra['bioproject'] == "PRJNA46333"][["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]


df_curated_PRJEB33701["Subject_ID"] = df_curated_PRJEB33701["Sample_accession"].apply(lambda x: sample_accession_subject_map[x])
df_curated_PRJEB33701["Nucleotide_Type"] = df_curated_PRJEB33701["Sample_accession"].apply(lambda x: sample_accession_nucleotide_map[x])
df_curated_PRJEB33701["Sex"] = df_curated_PRJEB33701["Sample_accession"].apply(lambda x: sample_accession_sex_map[x])
df_curated_PRJEB33701["Body_site"] = df_curated_PRJEB33701["Sample_accession"].apply(lambda x: sample_accession_body_site_map[x])
df_curated_PRJEB33701["Location"] = df_curated_PRJEB33701["Sample_accession"].apply(lambda x: sample_accession_location_map[x])
df_curated_PRJEB33701["Health_status"] = df_curated_PRJEB33701["Sample_accession"].apply(lambda x: sample_accession_health_map[x])

bins = [-1, 2, 5, 17, 39, 64, float('inf')]
labels = ['Infants', 'Child', 'Adolescent', 'Young_adult', 'Middle_age', 'Elderly']
df_curated_PRJEB33701["Age_catagories"] = pd.cut(
    df_curated_PRJEB33701["Age"], 
    bins=bins, 
    labels=labels, 
    right=True
)

df_curated_PRJEB33701["Age_catagories"] = df_curated_PRJEB33701["Age_catagories"].astype(object)

df_curated_PRJEB33701.to_csv("./Final_Curated_Out/PRJEB33701_Nguyen.csv", index = None)
