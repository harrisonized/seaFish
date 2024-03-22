The data is published in [this paper](https://doi.org/10.1038/s41598-020-70353-y):

```
Lantz C, Radmanesh B, Liu E, Thorp EB, Lin J. Single-cell RNA sequencing uncovers heterogenous transcriptional signatures in macrophages during efferocytosis. Sci Rep. 2020 Aug 31;10(1):14333. doi: 10.1038/s41598-020-70353-y. PMID: 32868786; PMCID: PMC7459098.
```
The data is publically available in GEO under the accession number [GSE156234](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156234).

Your directory should be structured the following way:

```
lantz-2020-efferocytosis/
├─ input/
│   ├─ WT_Ctrl/
│   │   └─ GSE156234_WT_Ctrl.txt.gz
│   ├─ WT_2_hr/
│   │   └─ GSE156234_WT_2_hr.txt.gz
│   ├─ WT_6_hr/
│   │   └─ GSE156234_WT_6_hr.txt.gz
│   ├─ MerKD_Ctrl/
│   │   └─ GSE156234_MerKD_Ctrl.txt.gz
│   ├─ MerKD_2hr/
│   │   └─ GSE156234_MerKD_2hr.txt.gz
│   └─ MerKD_6_hr/
│       └─ GSE156234_MerKD_6_hr.txt.gz
└─ README.md
```
In total, there are 6 files split from `GSE156234_aggregated_raw_counts.tsv.gz`.