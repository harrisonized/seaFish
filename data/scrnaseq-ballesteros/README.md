The data is published in [this paper](https://doi.org/10.1016/j.cell.2020.10.003):

```
2020, Ballesteros et al. Co-option of Neutrophil Fates by Tissue
Environments. https://doi.org/10.1016/j.cell.2020.10.003
```
The data is publically available in GEO under the accession number [GSE142754](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142754)

Your directory should be structured the following way:

```
scrnaseq-ballesteros/
├─ 10x_counts/
│   ├─ bm1/
│   │   └─ counts
│   │       ├─ GSM4239559_barcodes_BM1.tsv.gz
│   │       ├─ GSM4239559_BM1_mm10_genes.tsv.gz
│   │       └─ GSM4239559_BM1_mm10_matrix.mtx.gz
│   ├─ bm2/
│   │   └─ counts
│   │       ├─ GSM4239559_barcodes_BM1.tsv.gz
│   │       ├─ GSM4239559_BM1_mm10_genes.tsv.gz
│   │       └─ GSM4239559_BM1_mm10_matrix.mtx.gz
│   ├─ lung1/counts/...
│   ├─ lung2/counts/...
│   ├─ pbzt5/counts/...
│   ├─ pbzt13_1/counts/...
│   ├─ pbzt13_2/counts/...
│   ├─ spleen1/counts/...
│   └─ spleen2/counts/...
└─ README.md
```
In total, there are 27 files, split into 9 folders.