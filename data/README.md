# Pipeline input files

Excel exports ("Proteins" view with nested peptide rows) of the proteome-wide
Sequest HT searches against the UniProt human reference proteome, exported
from Proteome Discoverer 2.2. Folder and file names are kept exactly as
produced by the laboratory because `analyze_hspindly_ptms.py` refers to them.
Folder dates are export dates; acquisition dates are given below.

| Folder / file | Run | Search | Acquired |
|---------------|-----|--------|----------|
| `RPE-1 NZ y NZ+OK 22.02.22/NZ_vs_UniprotHuman-phospho.xlsx` | R1 | phospho | 2022-02-18 |
| `RPE-1 NZ y NZ+OK 22.02.22/NZ_vs_uniprotHuman-Ubi.xlsx` | R1 | GlyGly | 2022-02-18 |
| `RPE-1 NZ y NZ+OK 22.02.22/NZ-OK_vs_UniprotHuman-phospho.xlsx` | R2 | phospho | 2022-02-18 |
| `RPE-1 NZ y NZ+OK 22.02.22/NZ-OK_vs_UniprotHuman-Ubi.xlsx` | R2 | GlyGly | 2022-02-18 |
| `HCT116 NZ y NZ+OK 22.03.22/NZ_240322_vs_UniprotHuman-PHOSPHO.xlsx` | R3 | phospho | 2022-03-24 |
| `HCT116 NZ y NZ+OK 22.03.22/NZ_240322_vs_UniprotHuman-UBI.xlsx` | R3 | GlyGly | 2022-03-24 |
| `HCT116 NZ y NZ+OK 22.03.22/NZ-OK_240322_vs_UniprotHuman-PHOSPHO.xlsx` | R4 | phospho | 2022-03-24 |
| `HCT116 NZ y NZ+OK 22.03.22/NZ-OK_240322_vs_UniprotHuman-UBI.xlsx` | R4 | GlyGly | 2022-03-24 |
| `HCT116 NZ+PPasa y NZ+OK 08.02.22/NZ-OK-WashingBeadSol_vs_UniprotHuman_PTM-P.xlsx` | R5 | phospho | 2022-01-27 |
| `HCT116 NZ+PPasa y NZ+OK 08.02.22/NZ-OK-WashingBeadSol_vs_UniprotHuman_PTM-UBI.xlsx` | R5 | GlyGly | 2022-01-27 |
| `HCT116 NZ+PPasa y NZ+OK 08.02.22/NZPPasa-WashingBeadSol_vs_UniprotHuman_PTM-P.xlsx` | R6 | phospho | 2022-01-29 |
| `HCT116 NZ+PPasa y NZ+OK 08.02.22/NZPPasa-WashingBeadSol_vs_UniprotHuman_PTM-UBI.xlsx` | R6 | GlyGly | 2022-01-29 |

For R1, R2, R5 and R6 the same exports are also deposited in PRIDE PXD085140
under standardised names (SPDL1_R<n>_<date>_uniprot_<search>_<condition>.xlsx);
for R3 and R4 PRIDE holds the corresponding .pdResult files.
