The data is published in [this paper](https://doi.org/10.1016/j.celrep.2023.112064):

```
2023, Tomlinson et al. Staphylococcus aureus stimulates neutrophil
itaconate production that suppresses the oxidative burst.
https://doi.org/10.1016/j.celrep.2023.112064
```
The data is publically available in GEO under the accession number [GSE215195](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215195).

Your directory should be structured the following way:

```
tomlinson-2023-neutrophils-staph/
├─ input/
│   ├─ LAC-infected Irg1-null/
│   │   ├─ GSM6625635_AK4_barcodes.tsv.gz
│   │   ├─ GSM6625635_AK4_features.tsv.gz
│   │   ├─ GSM6625635_AK4_INDEX_SWAP_CONTAMINANT_READ_IDS.txt.gz
│   │   └─ GSM6625635_AK4_matrix.mtx.gz
│   ├─ LAC-infected WT/
│   │   ├─ GSM6625634_AK3_barcodes.tsv.gz
│   │   ├─ GSM6625634_AK3_features.tsv.gz
│   │   ├─ GSM6625634_AK3_INDEX_SWAP_CONTAMINANT_READ_IDS.txt.gz
│   │   └─ GSM6625634_AK3_matrix.mtx.gz
│   ├─ PBS-treated Irg1-null/...
│   └─ PBS-treated-WT/...
└─ README.md
```
In total, there are 16 files, split into 4 folders.