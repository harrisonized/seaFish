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

## Protocol

Animal studies. C57BL/6J mice were used as WT controls and bred in the Northwestern Center for Com- parative Medicine facility. MerKD (referred to herein as MerTK−/−) have been previously described25 and were backcrossed to C57BL/6J mice for ten generations. Eight- to twelve-week old mice were utilized for experiments. Mice were bred and housed in a pathogen-free, temperature- and humidity-controlled environment with access to standard mouse chow and water ad libitum.

**Ex vivo efferocytosis assay**

Resident peritoneal cells were harvested after lavage with cold saline. Peritoneal cells from 3 mice were pooled together for each experimental group. Initial cell selection was achieved through adherence to non-treated, low-adherence cell culture plates for 1 h and rinsed to remove non-adherent cells.

Apoptotic cells (ACs) were generated using GFP-Jurkat Human T cells (GenTarget) exposed to UV radiation for 7 min followed by a 2-h incubation at 37 °C. Apoptosis was confirmed by annexin V positive, propidium iodide negative identification affirming greater than 80% apoptosis.

Adherent resident peritoneal cells were co-cultivated with ACs at a ratio of 5 ACs to 1 peritoneal cell for 2 h or 6 h as indicated. Control cells were given a media change corresponding to the 6-h timepoint. Non-engulfed ACs were removed from the co-culture through rigorous rinsing with warm saline. Adherent macrophages were removed from the plate with Accutase and resuspended into a single cell suspension. Apoptotic cell engulfment was confirmed with confocal fluorescence microscopy.

**Single cell library preparation and RNA sequencing**

The single cell RNA-Seq libraries were prepared using the 10X Genomics Single Cell 5′ Gel Bead and Library Kit pipeline following manufacturer’s protocols. Cell suspensions were diluted to target a recovery of 4,000 cells per sample.
