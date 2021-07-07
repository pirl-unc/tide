This image is for running [TIDE](https://github.com/liulab-dfci/TIDEpy) for our
Nextflow modules.

## Decoding the tag structure
y.z 
y is the version of TIDE.  
z is the version of this Dockerfile.  

```bash
docker build . -t benjamin-vincent-lab/tide:1.1
docker run -it --rm benjamin-vincent-lab/tide:1.1 /bin/bash
```

```Python
import os
import pandas as pd
import numpy as np
from tidepy.pred import TIDE 

os.chdir('test')

df = pd.read_csv("entrez_rkpm_counts.tsv", sep='\t', header=None, index_col=False)

gene_names = df.iloc[0][1:].tolist()
gene_names = [int(i) for i in gene_names]
sample_names = df.iloc[1:, 0]
df = df.drop(columns=df.columns[0])
df = df.drop(labels=0, axis=0)

# transpose
df = df.transpose()
df = df.rename(index=lambda x: gene_names[x - 1])
df.columns = sample_names

# transform
df = df.add(1)
df = df.astype(float)
df = np.log2(df)

# now for each column subtract the column mean
df = df - df.mean()

tide_out_df = TIDE(df,cancer='Melanoma',pretreat=False)
tide_out_df.index.name = 'Run_ID'
tide_out_df.reset_index(inplace=True)

tide_out_df.columns = ('Run_ID', 'TIDE_No_Benefits', 'TIDE_Responder', 'TIDE', 'TIDE_IFNG', 'TIDE_MSI', 'TIDE_CD274','TIDE_CD8','TIDE_CTL_Flag','TIDE_Dysfunction','TIDE_Exclusion','TIDE_MDSC','TIDE_CAF','TIDE_TAM_M2','TIDE_CTL')

tide_out_df.to_csv('tide.tsv', sep = '\t', index = False)
```
