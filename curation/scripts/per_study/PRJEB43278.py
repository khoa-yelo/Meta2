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

df_curated_PRJEB43278 =  df_curated[df_curated['Study_ID'] == "PRJEB43278"]
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB43278"][["sample_accession", "sample_sample-name"]]
df_sra_subject = df_sra[df_sra['bioproject'] == "PRJDB5315"]
sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_nucleotide_map = df_sra_subject[["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]

df_curated_PRJEB43278["Subject_ID"] = df_curated_PRJEB43278["Sample_accession"].map(sample_accession_subject_map)
df_curated_PRJEB43278["Nucleotide_Type"] = df_curated_PRJEB43278["Sample_accession"].map(sample_accession_nucleotide_map)
df_curated_PRJEB43278["Body_site"] = "oral"
df_curated_PRJEB43278["Body_site_core"] = "oral"
df_curated_PRJEB43278["Nucleotide_Type"] = "DNA"

df_curated_PRJEB43278.to_csv("/content/drive/Shareddrives/Meta2/Metadata/Final_Curated_Out/PRJEB43278.csv", index = None)
