The data is published in [this paper](https://doi.org/10.1016/j.cell.2020.10.003):

```
2020, Ballesteros et al. Co-option of Neutrophil Fates by Tissue Environments.
https://doi.org/10.1016/j.cell.2020.10.003
```
The data is publically available in GEO under the accession number [GSE142754](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142754).

Your directory should be structured the following way:

```
ballesteros-2020/
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

## Protocol

Blood was taken through cardiac puncture with 1ml syringe attached to a 26 g needle filled with 50mL of 0.5M EDTA. Samples were lysed in RBC lysis buffer for 4min.

For BM cells, mice femurs were flushed using a 23–gauge needle in PBS containing 2mM EDTA and 2% fetal bovine serum (FBS) and passed through a 70–mm nylon mesh sieve.

Spleens were harvested and homogenized into single-cell suspensions using 70–mm nylon mesh sieves and syringe plungers.

Lungs were harvested and cut dry into small pieces before digested in liberase TM and DNase1 for 45min at 37C. Lungs were then homogenized into single-cell suspensions using 70–mm nylon mesh sieves and syringe plungers.

Blood, bone marrow, spleen and lung leukocyte suspensions from 6-12 week-old C57BL/6 male mice were FACS sorted in a BD Aria II sorter for CD45+CD11b+ cells. Cells were partitioned in Gel Beads in Emulsion (GEMs) and lysed, followed by RNA barcoding, reverse transcription and PCR amplification (12–14 cycles). Sequencing-ready scRNA-seq were prepared according to the manufacturer’s instructions, checked and quantified on the Agilent TapeStation 2200 and Qubit 3.0 instruments. Sequencing was performed on an Illumina NextSeq 500 machine using the NextSeq 500/550 High Output v2 kit (75 cycles), at the Center for Omics Sciences at the IRCCS Ospedale San Raffaele (COSR).
