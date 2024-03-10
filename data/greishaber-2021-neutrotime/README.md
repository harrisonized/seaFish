The data is published in [this paper](https://doi.org/10.1038/s41467-021-22973-9):

```
2021, Grieshaber-Bouyer R et al. The neutrotime transcriptional signature defines a single continuum of neutrophils across biological compartments.
https://doi.org/10.1038/s41467-021-22973-9
```
The data is publically available in GEO under the accession number [GSE165276](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165276).

Your directory should be structured the following way:

```
greishaber-2021-neutrotime/
├─ input/
│   ├─ blood-1/
│   │   └─ GSM5029335_BL_dataset1.txt.gz
│   ├─ blood-2/
│   │   └─ GSM5029338_BL_dataset2.txt.gz
│   ├─ blood-KBxN/
│   │   └─ GSM5029341_blood_KBxN.txt.gz
│   ├─ blood-unstim/
│   │   └─ GSM5029341_blood_unstim.txt.gz
│   ├─ bm-1/
│   │   └─ GSM5029336_BM_dataset1.txt.gz
│   ├─ bm-2/
│   │   └─ GSM5029339_BM_dataset2.txt.gz
│   ├─ joint-KBxN/
│   │   └─ GSM5029341_joint_KBxN.txt.gz
│   ├─ lung-IL1/
│   │   └─ GSM5029341_lung_IL1.txt.gz
│   ├─ peritoneum-IL1/
│   │   └─ GSM5029341_peritoneum_IL1.txt.gz
│   ├─ spleen-1/
│   │   └─ GSM5029337_SP_dataset1.txt.gz
│   └─ spleen-2/
│   │   └─ GSM5029340_SP_dataset2.txt.gz
└─ README.md
```
In total, there are 11 files, 6 are original and 5 are split from `GSM5029341_inflammation_dataset3.txt.gz` using 'tissue' and 'stimulus' metadata found in `GSM5029341_inflammation_dataset3_readme.txt.gz`.