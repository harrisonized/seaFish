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
│       └─ GSM5029340_SP_dataset2.txt.gz
└─ README.md
```
In total, there are 11 files, 6 are original and 5 are split from `GSM5029341_inflammation_dataset3.txt.gz` using 'tissue' and 'stimulus' metadata found in `GSM5029341_inflammation_dataset3_readme.txt.gz`.

## Protocol

Male WT C57BL/6 J mice were obtained from The Jackson Laboratory (Stock No 000664) at age of 5 weeks. Male and female K/BxN mice expressing the T-cell receptor transgene KRN and the major histocompatibility complex class II molecule Ag7 were housed in our animal facility at the Brigham and Women’s Hospital and serum from male and female mice was obtained at age of 8–11 weeks and pooled. Control mice were harvested at 6 weeks (dataset 1) and 8 weeks (dataset 2).

**Blood, Bone Marrow, and Spleen**

Neutrophils were obtained from the circulation of healthy anesthetized mice by cardiac puncture: 1 ml of blood was collected by cardiac puncture in a syringe coated with EDTA. Blood, bone marrow and spleen were obtained from the same mice in each experiment. Following cervical dislocation, mice were immediately dissected to obtain the spleen and the tibias + femurs. The spleen was carefully cleaned of any attached fat and lymph nodes and then minced in a cell culture dish with the sterile back of a syringe to dissociate splenic immune cells. The dissociated splenic tissue was then passed through a 70-micron filter to create a single cell suspension. Bone marrow from femurs and tibia was flushed using 4°C media to obtain bone marrow suspensions. All tissues were placed immediately into 4°C media.

**Experimental pneumonitis**

Male WT mice aged 8 weeks were anesthetized with 100 mg/kg ketamine and 16 mg/kg xylazine. 25 ng of recombinant Mouse IL-1β in 30 μL of phosphate-buffered saline were administrated intranasally. After 3 h, mice were euthanized and bronchoalveolar lavage (BAL) was performed using cold PBS.

**Experimental peritonitis**

Male WT mice aged 8 weeks were injected intraperitoneally with 25 ng IL-1β in 200 μL of PBS. After 3 h, mice were euthanized and peritoneal cells were harvested from the peritoneum with 5 mL of cold PBS.

**Neutrophil isolation**

All isolation solutions were at 4°C to avoid activation of neutrophils. 4°C phenol-red free DMEM + 0.1% sodium azide + 10 mM HEPES + 2% FCS + 5 mM EDTA media was used in all cell manipulation steps. Single-cell suspensions were gently pelleted at 400 × g for 5 min at 4°C and resuspended in 2 ml cold ACK Lysing Buffer for erythrocyte lysis. After 3 min of lysis, medium was added and cells were pelleted at 400 × g for 5 min at 4°C. Cells were then resuspended in cold media containing antibodies for staining.


**FACS**

Prior to sequencing, single-cell suspensions were stained with anti-mouse Ly6G (clone 1A8)-Alexa Fluor 647 and anti-mouse CD11b (clone M1/70)-Brilliant Violet 421TM at 1:100 dilution, for 30 min on ice. 10 min before fluorescence-activated cell sorting, propidium iodide was added to a final con- centration of 5 ng/ml. Neutrophils were gated based on FSC-A and SSC-A, doublets excluded in FSC-H vs. FSC-W and SSC-H vs. SSC-W, propidium iodide negative live cells and finally Ly6G-positive CD11b-positive neutrophils were selected. Cutoffs for sorting were determined based on unstained, isotype and single-stained cells. Cells were sorted directly into PBS with a final concentration of 0.04% BS. The sorting strategy is shown in Supplementary Fig. 1a. All steps were performed on ice with cold reagents, and total time from mouse euthanasia to single cell encapsulation was <2 h.

**Droplet-based single-cell RNA-sequencing.**

Sorted Ly6G-positive CD11b-positive neutrophils were loaded on a 10X Chromium device, following standard steps for library preparation, quality control, amplification and sequencing. 