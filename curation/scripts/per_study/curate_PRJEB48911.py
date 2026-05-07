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

df_curated_PRJEB48911 = df_curated[df_curated['Study_ID'] == "PRJEB48911"]
df_original_subject = df_original[df_original["study_bioproject"] == "PRJEB48911"][["sample_accession", "sample_sample-name","sample_sample-desc","sample_biosample"]]
df_sra_subject = df_sra[df_sra['bioproject'] == "PRJEB11419"]
sample_accession_subject_map = df_original_subject.set_index("sample_accession").to_dict()["sample_sample-name"]
sample_accession_nucleotide_map = df_sra_subject[["sample_accession", "nucleotide_type"]].set_index("sample_accession").to_dict()["nucleotide_type"]

df_sra[df_sra['bioproject'] == "PRJEB11419"] #use study_ID_original here
df_original_subject['sample_sample-desc'].unique()
df_original_subject = df_original_subject[df_original_subject['sample_sample-desc'] != 'American Gut control']

PRJEB48911_meta = pd.read_csv("/content/drive/Shareddrives/Meta2/Metadata/Manual_metadata/PRJEB48911.csv", sep = ',')
PRJEB48911_meta = PRJEB48911_meta.rename(columns={"BioSample" : 'sample_biosample'})
df_with_meta = df_original_subject.merge(PRJEB48911_meta, on="sample_biosample", how="left")
list(df_with_meta.columns)

mental_illness_map = {
    'Yes':         'Mental illness',
    'true':        'Mental illness',
    'No':          'No Mental illness',
    'false':       'No Mental illness',
    'Unspecified': np.nan,
}

df_with_meta['mental_illness'] = df_with_meta['mental_illness'].map(mental_illness_map)



antibiotic_history_map = {
    'Week':                                       'Antibiotics: Within past week',
    'Month':                                      'Antibiotics: Within past month',
    '6 months':                                   'Antibiotics: Within past 6 months',
    'Year':                                       'Antibiotics: Within past year',
    'I have not taken antibiotics in the past year.': 'Antibiotics: None in past year',
    'Unspecified':                                np.nan,
}

df_with_meta['antibiotic_history'] = df_with_meta['antibiotic_history'].map(antibiotic_history_map)
diagnosed_responses = {
    'Diagnosed by a medical professional (doctor, physician assistant)',
    'Self-diagnosed'
}

# Label map: column name → condition label
condition_labels = {
    'ibs': 'IBS',
    'ibd': 'IBD',
    'Diabetes': 'Diabetes',
    'autoimmune': 'Autoimmune',
    'cancer': 'Cancer',
    'cardiovascular_disease': 'Cardiovascular disease',
    'lung_disease': 'Lung disease',
    'kidney_disease': 'Kidney disease',
    'liver_disease': 'Liver disease',
    'thyroid': 'Thyroid disorder',
    'epilepsy_or_seizure_disorder': 'Epilepsy',
    'migraine': 'Migraine',
    'acid_reflux': 'Acid reflux',
    'skin_condition': 'Skin condition',
    'fungal_overgrowth': 'Fungal overgrowth',
    'cdiff': 'C. difficile',
    #'lactose': 'Lactose intolerance',
    #'gluten': 'Gluten sensitivity',
    'sibo': 'SIBO',
}

def recode_condition(value, label):
    if value in diagnosed_responses:
        return label
    elif value == 'I do not have this condition':
        return f'No {label}'
    else:  # nan or 'Unspecified'
        return np.nan

for col, label in condition_labels.items():
    df_with_meta[col] = df_with_meta[col].apply(lambda x: recode_condition(x, label))
df_with_meta['seasonal_allergies_map'].unique()
seasonal_allergies_map = {
    'Yes':         'lactose',
    'true':        'lactose',
    'No':          'No seasonal_allergies',
    'false':       'No seasonal_allergies',
    'Unspecified': np.nan,
}
df_with_meta['gluten'].unique()
df_with_meta['seasonal_allergies'] = df_with_meta['seasonal_allergies'].map(seasonal_allergies_map)

gluten_map = {
  
    'No':          'No gluten sensitivity',
    'Not provided':  np.nan,
    'Unspecified': np.nan,
}

df_with_meta['gluten'] = df_with_meta['gluten'].replace(gluten_map)

diet_cols = [
    'fruit_frequency', 'vegetable_frequency', 'whole_grain_frequency',
    'fermented_plant_frequency', 'probiotic_frequency', 'red_meat_frequency',
    'alcohol_frequency', 'sugar_sweetened_drink_frequency', 'artificial_sweeteners',
]

diet_labels = {
    'fruit_frequency': 'Fruit',
    'vegetable_frequency': 'Vegetables',
    'whole_grain_frequency': 'Whole grains',
    'fermented_plant_frequency': 'Fermented foods',
    'probiotic_frequency': 'Probiotics',
    'red_meat_frequency': 'Red meat',
    'alcohol_frequency': 'Alcohol',
    'sugar_sweetened_drink_frequency': 'Sugary drinks',
    'artificial_sweeteners': 'Artificial sweeteners',
    'exercise_frequency' : 'exercise', 
    'smoking_frequency' : 'smoking'
}

frequency_map = {
    'Daily':                        lambda l: f'Daily {l}',
    'Regularly (3-5 times/week)':   lambda l: f'Regular {l}',
    'Occasionally (1-2 times/week)': lambda l: f'Occasional {l}',
    'Rarely (a few times/month)':   lambda l: f'Rare {l}',
    'Never':                        lambda l: f'No {l}',
    'Unspecified':                  lambda l: np.nan,
}

def recode_frequency(value, label):
    if value in frequency_map:
        return frequency_map[value](label)
    return np.nan  # catches nan and any unexpected values

for col, label in diet_labels.items():
    df_with_meta[col] = df_with_meta[col].apply(lambda x: recode_frequency(x, label))

cancer_treatment_map = {
    'Surgery only':       'Cancer treatment: Surgery only',
    'Chemotherapy':       'Cancer treatment: Chemotherapy',
    'Radiation therapy':  'Cancer treatment: Radiation therapy',
    'No treatment':       'Cancer treatment: None',
    'Unspecified':        np.nan,
    'Not provided':       np.nan,
}

df_with_meta['cancer_treatment'] = df_with_meta['cancer_treatment'].map(cancer_treatment_map)

sleep_map = {
    'Less than 5 hours': 'Sleep: Less than 5 hours',
    '5-6 hours':         'Sleep: 5-6 hours',
    '6-7 hours':         'Sleep: 6-7 hours',
    '7-8 hours':         'Sleep: 7-8 hours',
    '8 or more hours':   'Sleep: 8+ hours',
    'Unspecified':       np.nan,
}

df_with_meta['sleep_duration'] = df_with_meta['sleep_duration'].map(sleep_map)

animal_treatment_map = {
    'Yes':       'consume animal product',
    'No':       'doesn\'t consume animal product',
    'true':  'consume animal product',
    'Unspecified':        np.nan,
    'Not sure':       np.nan,
    'false' : 'doesn\'t consume animal product'
}

df_with_meta['consume_animal_products_abx'] = df_with_meta['consume_animal_products_abx'].map(animal_treatment_map)

df_with_meta['sample_sample-desc'].unique()
body_site_map = {
    'American Gut Project Stool sample':        ('Gut',      'Gut'),
    'American Gut Project Stool Sample':        ('Gut',      'Gut'),  # capital S variant
    'American Gut Project Mouth sample':        ('Oral',     'Oral'),
    'American Gut Project Nasal mucus sample':  ('Nasal',    'Nasal'),
    'American Gut Project Nares sample':        ('Nares',    'Nasal'),
    'American Gut Project Left Hand sample':    ('Hand',     'Skin'),
    'American Gut Project Forehead sample':     ('Forehead', 'Skin'),
    'American Gut control':                      (np.nan,    np.nan),
}

df_with_meta['body_site']      = df_with_meta['sample_sample-desc'].map(lambda x: body_site_map.get(x, (np.nan, np.nan))[0])
df_with_meta['body_site_core'] = df_with_meta['sample_sample-desc'].map(lambda x: body_site_map.get(x, (np.nan, np.nan))[1])

health_cols = [
    'ibs', 'ibd', 'Diabetes', 'autoimmune', 'cancer',
    'cardiovascular_disease', 'lung_disease', 'kidney_disease',
    'liver_disease', 'thyroid', 'mental_illness',
    'epilepsy_or_seizure_disorder', 'migraine', 'acid_reflux',
    'skin_condition', 'fungal_overgrowth', 'cdiff', 'seasonal_allergies',
    'gluten', 'sibo', 'antibiotic_history', 'cancer_treatment','ibd_diagnosis_refined'
] #left lactose out 

# Combine all conditions into one column (already recoded columns from earlier)
def combine_health(row):
    conditions = [row[col] for col in health_cols if pd.notna(row[col])]
    return ', '.join(conditions) if conditions else np.nan

df_with_meta['Health_status'] = df_with_meta.apply(combine_health, axis=1)


diet_cols = ['diet_type', 'alcohol_frequency','fruit_frequency', 'vegetable_frequency', 'whole_grain_frequency',
    'fermented_plant_frequency', 'probiotic_frequency', 'red_meat_frequency',
    'sugar_sweetened_drink_frequency', 'artificial_sweeteners',
    'consume_animal_products_abx', 'gluten']
def combine_diet(row):
  conditions = [row[col] for col in diet_cols if pd.notna(row[col])]
  return ', '.join(conditions) if conditions else np.nan

df_with_meta['diet'] = df_with_meta.apply(combine_diet, axis=1)


lifestyle_cols = ['exercise_frequency', 'smoking_frequency', 'sleep_duration']
def combine_lifestyle(row):
  conditions = [row[col] for col in lifestyle_cols if pd.notna(row[col])]
  return ', '.join(conditions) if conditions else np.nan

df_with_meta['lifestyle'] = df_with_meta.apply(combine_lifestyle, axis=1)

df_with_meta.rename(columns={'sample_accession': 'Sample_accession'}, inplace=True)
  
sample_accession_sex_map = (df_with_meta[["Sample_accession", "sex"]].set_index("Sample_accession").to_dict()["sex"])
sample_accession_health_map = (df_with_meta[["Sample_accession", "Health_status"]].set_index("Sample_accession").to_dict()["Health_status"])
sample_accession_location_map = (df_with_meta[["Sample_accession", "country_of_birth"]].set_index("Sample_accession").to_dict()["country_of_birth"])
sample_accession_agecat_map= (df_with_meta[["Sample_accession", "age_cat"]].set_index("Sample_accession").to_dict()["age_cat"])
sample_accession_age_map = (df_with_meta[["Sample_accession", "Age"]].set_index("Sample_accession").to_dict()["Age"])
sample_accession_diet_map = (df_with_meta[["Sample_accession", "diet"]].set_index("Sample_accession").to_dict()["diet"])
sample_accession_lifestyle_map = (df_with_meta[["Sample_accession", "lifestyle"]].set_index("Sample_accession").to_dict()["lifestyle"])
sample_accession_body_site_map = (df_with_meta[["Sample_accession", "body_site"]].set_index("Sample_accession").to_dict()["body_site"])
sample_accession_body_site_core = (df_with_meta[["Sample_accession", "body_site_core"]].set_index("Sample_accession").to_dict()["body_site_core"])

df_curated_PRJEB48911["Subject_ID"] = df_curated_PRJEB48911["Sample_accession"].map(sample_accession_subject_map)
df_curated_PRJEB48911["Nucleotide_Type"] = 'DNA'
df_curated_PRJEB48911["Body_site"] = df_curated_PRJEB48911["Sample_accession"].map(sample_accession_body_site_map)
df_curated_PRJEB48911["Body_site_core"] = df_curated_PRJEB48911["Sample_accession"].map(sample_accession_body_site_core)
df_curated_PRJEB48911["Age"] = df_curated_PRJEB48911["Sample_accession"].map(sample_accession_age_map)
df_curated_PRJEB48911["Age_cat"] = df_curated_PRJEB48911["Sample_accession"].map(sample_accession_agecat_map)
df_curated_PRJEB48911["Sex"] = df_curated_PRJEB48911["Sample_accession"].map(sample_accession_sex_map)
df_curated_PRJEB48911["Health_status"] = df_curated_PRJEB48911["Sample_accession"].map(sample_accession_health_map)
df_curated_PRJEB48911["Lifestyle"] = df_curated_PRJEB48911['Sample_accession'].map(sample_accession_lifestyle_map)
df_curated_PRJEB48911["Diet"] = df_curated_PRJEB48911['Sample_accession'].map(sample_accession_diet_map)
df_curated_PRJEB48911["Location"] = df_curated_PRJEB48911['Sample_accession'].map(sample_accession_location_map)

df_curated_PRJEB48911

df_curated_PRJEB48911.to_csv("/content/drive/Shareddrives/Meta2/Metadata/Final_Curated_Out/PRJEB48911.csv", index = None)
