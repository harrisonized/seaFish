The data is published in [this paper](https://doi.org/10.1016/j.cell.2020.10.003):

```
2020, Ballesteros et al. Co-option of Neutrophil Fates by Tissue
Environments. https://doi.org/10.1016/j.cell.2020.10.003
```
The data is publically available in GEO under the accession number [GSE142754](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142754)

Your directory should be structured the following way:

```
scrnaseq-ballesteros/
├─ input/
│   ├─ bm1/
│   │   ├─ GSM4239559_barcodes_BM1.tsv.gz
│   │   ├─ GSM4239559_BM1_mm10_genes.tsv.gz
│   │   └─ GSM4239559_BM1_mm10_matrix.mtx.gz
│   ├─ bm2/
│   │   ├─ GSM4239559_barcodes_BM1.tsv.gz
│   │   ├─ GSM4239559_BM1_mm10_genes.tsv.gz
│   │   └─ GSM4239559_BM1_mm10_matrix.mtx.gz
│   ├─ lung1/...
│   ├─ lung2/...
│   ├─ pbzt5/...
│   ├─ pbzt13_1/...
│   ├─ pbzt13_2/...
│   ├─ spleen1/...
│   ├─ spleen1/...
│   └─ config.json
└─ README.md
```
In total, there are 27 files, split into 9 folders. Note that there's a bug in `pbzt13_2` that prevents it from being used. The matrix contains 3,388 genes, but the genes file only contains 3,204 genes.