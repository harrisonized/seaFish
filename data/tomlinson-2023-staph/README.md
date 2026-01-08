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

## Protocol

Breeding pairs of Irg1-/- C57BL/6NJ mice were purchased from Jackson Laboratories and used to establish a colony at Columbia University. Wildtype (WT) C57BL/6NJ mice were purchased from Jackson Laboratories for each experiment.

8-week old WT and Irg1-/- C57BL/6NJ mice were intranasally infected with 2 x 10^7 USA300 LAC in 50 mL PBS or treated with PBS alone. 24 h after infection, mice were euthanized and their lungs were harvested. To create a single cell suspension, the lungs were placed in an Eppendorf tube containing an enzymatic digestion solution of 2 mg/mL collagenase I, 20 mg/mL dispase, 1 mg/mL elastase, and 1 mL/mL DNAse in PBS. The lungs were minced in the Eppendorf tube, then incubated with shaking at 37°C for 30 min. The digestion was quenched with 4 volumes of PBS supplemented with 10% hiFBS, and strained over a 70 mm cell strainer. The cell suspension was centrifuged at 1400 x g for 7 min at 4°C. Red blood cell lysis was performed with the eBioscience RBC lysis buffer.

The resulting cell pellet was resuspended in PBS supplemented with 0.04% bovine serum albumin before being loaded onto the 10X Genomics Chromium Single Cell Controller. Cell viability was above 95% for each sample. Approximately 5000 cells were analyzed per sample.

