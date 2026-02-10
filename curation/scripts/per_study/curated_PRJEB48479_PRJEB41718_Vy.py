from google.colab import drive    
drive.mount('/content/drive')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sra_cleaned_path = "/content/drive/Shareddrives/Meta2/Metadata/cleaned_sra_metadata.tsv"  #insert file used for curated
curated_cleaned_path = "/content/drive/Shareddrives/Meta2/Metadata/curated_metadata.csv"   
original_metadata_path = "/content/drive/Shareddrives/Meta2/Data/hmb_assemblies_metadata.csv"   
df_sra = pd.read_csv(sra_cleaned_path, sep="\t")
df_curated = pd.read_csv(curated_cleaned_path, sep=",")
df_original = pd.read_csv(original_metadata_path, sep=",")
df_curated.head()  #check the first 5 rows
print(df_curated.columns)   #print out columns of this file
df_curated[df_curated["Study_ID"] == "PRJEB48479"]
df_original.head()
df_original[df_original["study_bioproject"] == "PRJEB48479"]["sample_species"].value_counts(dropna=False)
df_curated[df_curated["Study_ID"] == "PRJEB48479"]["Body_site"].value_counts(dropna=False)
df_curated[df_curated["Study_ID"] == "PRJEB48479"]["Sex"].value_counts(dropna=False)
df_curated[df_curated["Study_ID"] == "PRJEB48479"]["Body_site_core"].value_counts(dropna=False)
df_curated[df_curated["Study_ID"] == "PRJEB48479"]["Health_status"].value_counts(dropna=False)
df_curated[df_curated["Study_ID"] == "PRJEB48479"]["Diet"].value_counts(dropna=False)
df_curated[df_curated["Study_ID"] == "PRJEB48479"]["Location"].value_counts(dropna=False)
df_curated[df_curated["Study_ID"] == "PRJEB48479"]["Lifestyle"].value_counts(dropna=False)
df_curated[df_curated["Study_ID"] == "PRJEB48479"]["Age"].value_counts(dropna=False)
for col in df_original.columns:
  print(col)
  
df_original[df_original["study_bioproject"] == "PRJEB48479"]["sample_geo-loc-name"].value_counts(dropna=False)
df_original[df_original["study_bioproject"] == "PRJEB48479"]["sample_samplemeta_body_mass_index"].value_counts(dropna=False)
df_original[df_original["study_bioproject"] == "PRJEB48479"]["sample_samplemeta_age"].value_counts(dropna=False)
for col in df_sra.columns:
  print (col)  
  # Used for each Study_ID
