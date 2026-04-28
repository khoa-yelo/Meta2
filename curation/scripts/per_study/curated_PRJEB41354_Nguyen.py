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


df_curated_PRJEB41354 = df_curated[df_curated['Study_ID'] == "PRJEB41354"].copy()

#create file
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB41354"][["sample_accession", "sample_sample-name"]]
df_original_subject["sample_sample-name"] = df_original_subject["sample_sample-name"].apply(lambda x: x)

sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_body_site_map = df_sra[df_sra['bioproject'] == "PRJEB28422"][["sample_accession", "body_site"]].set_index("sample_accession").to_dict()["body_site"]
sample_accession_age_map = df_sra[df_sra['bioproject'] == "PRJEB28422"][["sample_accession", "age_years"]].set_index("sample_accession").to_dict()["age_years"]
sample_accession_location_map = df_sra[df_sra['bioproject'] == "PRJEB28422"][["sample_accession", "geo_location"]].set_index("sample_accession").to_dict()["geo_location"]
sample_accession_nucleotide_map = df_sra[df_sra['bioproject'] == "PRJEB28422"][["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]

df_curated_PRJEB41354["Subject_ID"] = df_curated_PRJEB41354["Sample_accession"].apply(lambda x: sample_accession_subject_map[x])
df_curated_PRJEB41354["Nucleotide_Type"] = df_curated_PRJEB41354["Sample_accession"].apply(lambda x: sample_accession_nucleotide_map[x])
df_curated_PRJEB41354["Body_site"] = df_curated_PRJEB41354["Sample_accession"].apply(lambda x: sample_accession_body_site_map[x])
df_curated_PRJEB41354["Body_site_core"] = "Oral"
df_curated_PRJEB41354["Age"] = df_curated_PRJEB41354["Sample_accession"].apply(lambda x: sample_accession_age_map[x])
df_curated_PRJEB41354["Location"] = df_curated_PRJEB41354["Sample_accession"].apply(lambda x: sample_accession_location_map[x])

#supp
supp_path = "PRJEB41354_supp_Nguyen.xlsx"
df_supp = pd.read_excel(supp_path)
df_supp.columns.values[0] = "ID"
df_supp["ID"] = df_supp["ID"].astype(str)
df_supp["replace_specI"] = df_supp["replace_specI"].astype(str)
df_curated_PRJEB41354["Subject_ID"] = df_curated_PRJEB41354["Subject_ID"].astype(str)

sex_id_map = df_supp.set_index("ID")["sex"].to_dict()
group_id_map = df_supp.set_index("ID")["group"].to_dict()
sex_replace_map = df_supp.set_index("replace_specI")["sex"].to_dict()
group_replace_map = df_supp.set_index("replace_specI")["group"].to_dict()
age_replace_map = df_supp.set_index("replace_specI")["age"].to_dict()

sex_val = df_curated_PRJEB41354["Subject_ID"].map(sex_id_map)
group_val = df_curated_PRJEB41354["Subject_ID"].map(group_id_map)
df_curated_PRJEB41354["Sex"] = sex_val.fillna(df_curated_PRJEB41354["Subject_ID"].map(sex_replace_map))
df_curated_PRJEB41354["Health_status"] = group_val.fillna(df_curated_PRJEB41354["Subject_ID"].map(group_replace_map))

if "Age" not in df_curated_PRJEB41354.columns:
    df_curated_PRJEB41354["Age"] = np.nan

mask_age_missing = df_curated_PRJEB41354["Age"].isna()
df_curated_PRJEB41354.loc[mask_age_missing, "Age"] = df_curated_PRJEB41354.loc[mask_age_missing, "Subject_ID"].map(age_replace_map)
df_curated_PRJEB41354["Sex"] = df_curated_PRJEB41354["Sex"].replace({"F": "Female", "M": "Male"})
df_curated_PRJEB41354["Age"] = pd.to_numeric(df_curated_PRJEB41354["Age"], errors='coerce')


mask_germany_missing_health = (df_curated_PRJEB41354["Location"] == "Germany") & (df_curated_PRJEB41354["Health_status"].isna())
df_curated_PRJEB41354.loc[mask_germany_missing_health, "Health_status"] = "Control"
#age_catagories
bins = [-1, 2, 5, 17, 39, 64, float('inf')]
labels = ['Infants', 'Child', 'Adolescent', 'Young_adult', 'Middle_age', 'Elderly']

df_curated_PRJEB41354["Age_catagories"] = pd.cut(
    df_curated_PRJEB41354["Age"], 
    bins=bins, 
    labels=labels, 
    right=True
)

df_curated_PRJEB41354["Age_catagories"] = df_curated_PRJEB41354["Age_catagories"].astype(object)
df_curated_PRJEB41354.to_csv("./Final_Curated_Out/PRJEB41354_Nguyen.csv", index = None)
