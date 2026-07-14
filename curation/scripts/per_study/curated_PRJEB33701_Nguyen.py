# mount drive
from google.colab import drive
drive.mount('/content/drive')


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)

sra_cleaned_path = "/content/drive/Shareddrives/Meta2/Metadata/Curated_org_data/cleaned_sra_metadata.tsv"
curated_cleaned_path = "/content/drive/Shareddrives/Meta2/Metadata/Curated_org_data/curated_metadata.csv"
original_metadata_path = "/content/drive/Shareddrives/Meta2/Data/hmb_assemblies_metadata.csv"

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
df_curated_PRJEB33701["Health_status"] = "Healthy"


df_curated_PRJEB33701["Age_catagories"] = "Young_adult"

#body_site

supp_path = "PRJEB33701_supp_Nguyen.xlsx"
df_supp_s1 = pd.read_excel(supp_path, sheet_name="S1", skiprows=3, engine="openpyxl")
df_supp_s1.columns = df_supp_s1.columns.astype(str).str.strip()
df_supp_s1["SAMPLE"] = df_supp_s1["SAMPLE"].astype(str)
df_curated_PRJEB33701["Subject_ID"] = df_curated_PRJEB33701["Subject_ID"].astype(str)
sample_site_map = df_supp_s1.set_index("SAMPLE")["Site-Symmetry"].to_dict()
df_curated_PRJEB33701["raw_site"] = df_curated_PRJEB33701["Subject_ID"].map(sample_site_map)
body_site_dict = {
    "Gb": "Glabella", "Ea": "External auditory canal", "Na": "Nare",
    "Mb": "Manubrium", "Ac": "Antecubital fossa", "Vf": "Volar forearm",
    "Hp": "Hypothenar palm", "Ic": "Inguinal crease", "Tw": "Toe web space",
    "Ch": "Cheek", "Al": "Alar crease", "Ra": "Retroauricular crease",
    "Oc": "Occiput", "Ba": "Back", "Id": "Interdigital web",
    "Pc": "Popliteal fossa", "Tn": "Toenail", "Ph": "Plantar heel"
}
side_dict = {"R": "Right", "L": "Left"}
split_site = df_curated_PRJEB33701["raw_site"].str.split("-", expand=True)
site_eng = split_site[0].map(body_site_dict).fillna(split_site[0])

if 1 in split_site.columns:
    clean_side = split_site[1].str.split(":", expand=True)[0]
    side_eng = clean_side.map(side_dict)
else:
    
    side_eng = pd.Series(index=df_curated_PRJEB33701.index, dtype=object)

df_curated_PRJEB33701["Body_site"] = site_eng
mask_has_side = side_eng.notna()
df_curated_PRJEB33701.loc[mask_has_side, "Body_site"] = side_eng[mask_has_side] + " " + site_eng[mask_has_side]
df_curated_PRJEB33701 = df_curated_PRJEB33701.drop(columns=["raw_site"])

#sex
supp_path = "PRJEB33701_supp_Nguyen.xlsx"
df_supp_s1 = pd.read_excel(supp_path, sheet_name="S1", skiprows=3, engine="openpyxl")
df_supp_s1.columns = df_supp_s1.columns.astype(str).str.strip()
df_supp_s1["SAMPLE"] = df_supp_s1["SAMPLE"].astype(str)
df_curated_PRJEB33701["Subject_ID"] = df_curated_PRJEB33701["Subject_ID"].astype(str)
sample_gender_map = df_supp_s1.set_index("SAMPLE")["Gender"].to_dict()
df_curated_PRJEB33701["Sex"] = df_curated_PRJEB33701["Subject_ID"].map(sample_gender_map)
df_curated_PRJEB33701["Sex"] = df_curated_PRJEB33701["Sex"].replace({
    "F": "Female", 
    "M": "Male", 
    "female": "Female", 
    "male": "Male"
})


df_curated_PRJEB33701.to_csv("./Final_Curated_Out/PRJEB33701_Nguyen.csv", index = None)
