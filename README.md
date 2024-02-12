This image is for running [TIDE](https://github.com/liulab-dfci/TIDEpy) for our
Nextflow modules.

* this repository originally resided at https://github.com/Benjamin-Vincent-Lab/tide , it was moved here on 2/12/2024

## Decoding the tag structure
y.z 
y is the version of TIDE.  
z is the version of this Dockerfile.  

```bash
docker build . -t benjamin-vincent-lab/tide:1.2
docker run -it --rm benjamin-vincent-lab/tide:1.2 /bin/bash
```

```Python
#!/usr/bin/env python3
  
import pandas as pd
import numpy as np
from tidepy.pred import TIDE 
import glob

entrez_rkpm_path = "${entrez_rkpm}"
sample_data_path = "${sample_data}"

key_col = "Run_ID"
dataset_col = "Dataset"
tcga_study_col = "TCGA_Study"
prior_tx_col = "Prior_ICI_Tx"

df = pd.read_csv(entrez_rkpm_path, sep='\t', header=None, index_col=False)
sample_df = pd.read_csv(sample_data_path, sep='\t', usecols=[key_col,dataset_col,tcga_study_col,prior_tx_col])

gene_names = df.iloc[0][1:].tolist()
gene_names = [int(i) for i in gene_names]
sample_names = df.iloc[1:, 0]
df = df.drop(columns=df.columns[0])
df = df.drop(labels=0, axis=0)
df = df.transpose()
df = df.rename(index=lambda x: gene_names[x - 1])
df.columns = sample_names

df = df.astype(float)

combined_df = pd.DataFrame()
datasets = sample_df[dataset_col].unique()
for dataset in datasets:
  dataset_df = sample_df.loc[sample_df[dataset_col] == dataset]
  study_types = dataset_df[tcga_study_col].unique()
  for study_type in study_types:
    if study_type == "SKCM":
      tide_cancer_type = "Melanoma"
    elif study_type == "LUAD" or study_type == "LUSC":
      # 'both LUAD and LUSC, belong to the family of NSCLCs': https://www.nature.com/articles/s41598-020-77284-8
      tide_cancer_type = "NSCLC"
    else:
      tide_cancer_type = "Other"
    
    tissue_df = dataset_df.loc[dataset_df[tcga_study_col] == study_type]
    pretreated_ids = tissue_df.loc[tissue_df[prior_tx_col] == True][key_col].tolist()
    pretreated_ids = list(set(pretreated_ids).intersection(set(df.columns)))
    
    non_pretreated_ids = tissue_df.loc[tissue_df[prior_tx_col] != True][key_col].tolist()
    non_pretreated_ids = list(set(non_pretreated_ids).intersection(set(df.columns)))
    
    if len(pretreated_ids) > 0:
      tide_out_df = TIDE(df[pretreated_ids],cancer=tide_cancer_type, pretreat=True, ignore_norm=False, force_normalize=True)
      combined_df = combined_df.append(tide_out_df)
      del tide_out_df
    if len(non_pretreated_ids) > 0:
      tide_out_df = TIDE(df[non_pretreated_ids],cancer=tide_cancer_type, pretreat=False, ignore_norm=False, force_normalize=True)
      combined_df = combined_df.append(tide_out_df)
      del tide_out_df
              
combined_df.index.name = 'Run_ID'
combined_df.reset_index(inplace=True)
combined_df.columns = ('Run_ID', 'TIDE_No_Benefits', 'TIDE_Responder', 'TIDE', 'TIDE_IFNG', 'TIDE_MSI','TIDE_CD274','TIDE_CD8','TIDE_CTL_Flag','TIDE_Dysfunction','TIDE_Exclusion','TIDE_MDSC','TIDE_CAF','TIDE_TAM_M2','TIDE_CTL')
  
combined_df.to_csv('tide.tsv', sep = '\t', index = False)
```
