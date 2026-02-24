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
df_curated_PRJEB41718 =  df_curated[df_curated['Study_ID'] == "PRJEB41718"]
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB41718"][["sample_accession", "sample_sample-name"]]
def extract_subject_id(x):
    if pd.isna(x):      #fix lỗi split"_" nếu trong giá trị không có "_"
        return np.nan

    parts = str(x).split("_")
    if len(parts) > 1:
        return parts[1]
    else:
        return np.nan
df_original_subject["sample_sample-name"] = df_original_subject["sample_sample-name"].apply(extract_subject_id)
sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_nucleotide_map = df_sra[df_sra['bioproject'] == "PRJEB6070"][["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]
sample_accession_age_map = df_sra[df_sra['bioproject'] == "PRJEB6070"][["sample_accession","age_years"]].set_index("sample_accession").to_dict()["age_years"]  #dùng sample_accession làm key để map thông tin
df_curated_PRJEB41718["Subject_ID"] = df_curated_PRJEB41718["Sample_accession"].apply(lambda x: sample_accession_subject_map[x])
df_curated_PRJEB41718["Nucleotide_Type"] = df_curated_PRJEB41718["Sample_accession"].apply(lambda x: sample_accession_nucleotide_map[x])  #cột được thêm vào
df_curated_PRJEB41718["Age"] = df_curated_PRJEB41718["Sample_accession"].apply(lambda x: sample_accession_age_map[x])     #cột được thêm vào
df_curated_PRJEB41718["Location"] = "France"
df_curated_PRJEB41718.to_csv("/content/drive/Shareddrives/Meta2/Metadata/Final_Curated_Out/PRJEB41718.csv", index = None)