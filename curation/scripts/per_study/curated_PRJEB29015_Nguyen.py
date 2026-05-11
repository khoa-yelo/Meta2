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


df_curated_PRJEB29015 = df_curated[df_curated['Study_ID'] == "PRJEB29015"].copy()

#create file
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB29015"][["sample_accession", "sample_sample-name"]]
df_original_subject["sample_sample-name"] = df_original_subject["sample_sample-name"].apply(lambda x: x)

sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_location_map = df_sra[df_sra['bioproject'] == "PRJEB6997"][["sample_accession", "geo_location"]].set_index("sample_accession").to_dict()["geo_location"]
sample_accession_nucleotide_map = df_sra[df_sra['bioproject'] == "PRJEB6997"][["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]
sample_accession_body_site_map = df_original[df_original['study_bioproject'] == "PRJEB29015"].set_index("sample_accession")["sample_sample-desc"].to_dict()


df_curated_PRJEB29015["Subject_ID"] = df_curated_PRJEB29015["Sample_accession"].apply(lambda x: sample_accession_subject_map[x])
df_curated_PRJEB29015["Nucleotide_Type"] = df_curated_PRJEB29015["Sample_accession"].apply(lambda x: sample_accession_nucleotide_map[x])
df_curated_PRJEB29015["Body_site"] = df_curated_PRJEB29015["Sample_accession"].map(sample_accession_body_site_map)
df_curated_PRJEB29015["Body_site_core"] = np.nan
df_curated_PRJEB29015["Body_site_core"] = df_curated_PRJEB29015["Body_site_core"].astype(object)
df_curated_PRJEB29015.loc[df_curated_PRJEB29015["Body_site"].str.contains("oral", case=False, na=False), "Body_site_core"] = "oral"
df_curated_PRJEB29015.loc[df_curated_PRJEB29015["Body_site"].str.contains("Stool", case=False, na=False), "Body_site_core"] = "gut"
df_curated_PRJEB29015["Age_catagories"] = np.nan
df_curated_PRJEB29015["Age_catagories"] = df_curated_PRJEB29015["Age_catagories"].astype(object)


#sex_supp
supp_path = "PRJEB29015_supp_Nguyen.xlsx"
df_supp = pd.read_excel(supp_path, skiprows=2)
df_supp["NCBI BioSample ID"] = df_supp["NCBI BioSample ID"].astype(str)

biosample_sex_map = df_supp.set_index("NCBI BioSample ID")["Sex (F/M)"].to_dict()
df_original["sample_biosample"] = df_original["sample_biosample"].astype(str)
sample_accession_biosample_map = df_original[df_original['study_bioproject'] == "PRJEB29015"].set_index("sample_accession")["sample_biosample"].to_dict()
biosample_temp = df_curated_PRJEB29015["Sample_accession"].map(sample_accession_biosample_map)
df_curated_PRJEB29015["Sex"] = biosample_temp.map(biosample_sex_map)
df_curated_PRJEB29015["Sex"] = df_curated_PRJEB29015["Sex"].replace({
    "F": "Female", 
    "M": "Male"
})

#health_sup

biosample_health_map = df_supp.set_index("NCBI BioSample ID")["Subject type"].to_dict()
df_original["sample_biosample"] = df_original["sample_biosample"].astype(str)
sample_accession_biosample_map = df_original[df_original['study_bioproject'] == "PRJEB29015"].set_index("sample_accession")["sample_biosample"].to_dict()
biosample_temp = df_curated_PRJEB29015["Sample_accession"].map(sample_accession_biosample_map)
df_curated_PRJEB29015["Health_status"] = biosample_temp.map(biosample_health_map)
health_status_rename = {
    "Treated": "Treated Rheumatoid Arthritis patients",
    "Untreated": "Untreated Rheumatoid Arthritis patients"
}
df_curated_PRJEB29015["Health_status"] = df_curated_PRJEB29015["Health_status"].replace(health_status_rename)

df_curated_PRJEB29015.to_csv("./Final_Curated_Out/PRJEB29015_Nguyen.csv", index = None)
