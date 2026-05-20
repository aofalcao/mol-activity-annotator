# CSANNO readable report

## Metadata

| Field | Value |
| --- | --- |
| Query date/time | 2026-05-20T16:56:31+01:00 |
| Input SMILES | CN(C)C(=O)Cc1c(nc2ccc(C)cn12)c3ccc(C)cc3 |
| Molecule count | 1 |
| CSANNO version | CSANNO - (C) 2019/2026 - Andre O. Falcao BioISI - DI/FCUL version 0.5.20260515 |
| ChEMBL server | https://www.ebi.ac.uk |
| ChEMBL request warnings | none |

## CSANNO options

**Command used:**

```bash
csanno.py -mol 'CN(C)C(=O)Cc1c(nc2ccc(C)cn12)c3ccc(C)cc3' -sim 0.7 -out 'data\zolpidem-A'
```

- Search evidence included: active / positive evidence.
- Report sections requested: aggregate target-frequency tables; molecule-by-molecule detail tables.
- Target labels shown as: gene symbol.
- Human-readable report format: markdown.
- Tier 1 enabled: similar ChEMBL molecules were retrieved and self matches were excluded; the similarity threshold was 0.700.
- Activity classification used rule profile 'csanno-default-0.5' from C:\Users\aofal\Desktop\Projectos\CSANNO\csanno_default_rules.yaml.

## Raw data

### Tier 0 aggregate - exact ChEMBL matches

| Target | Molecule count | Fraction |
| --- | --- | --- |
| ADRA1A | 1 | 1.0000 |
| CHRM1 | 1 | 1.0000 |
| GABRA1 | 1 | 1.0000 |
| GABRA2 | 1 | 1.0000 |
| GABRA3 | 1 | 1.0000 |
| HDAC6 | 1 | 1.0000 |
| LMNA | 1 | 1.0000 |
| REP | 1 | 1.0000 |
| SLCO1B1 | 1 | 1.0000 |
| SLCO1B3 | 1 | 1.0000 |
| TSPO | 1 | 1.0000 |

### Tier 0 detail - exact ChEMBL matches by input molecule

| Molecule | Annotation count | Annotations |
| --- | --- | --- |
| User1 | 11 | ADRA1A, CHRM1, GABRA1, GABRA2, GABRA3, HDAC6, LMNA, REP, SLCO1B1, SLCO1B3, TSPO |

### Tier 1 aggregate — similar ChEMBL analogue pool

No annotations found.

### Tier 1 detail — analogue-derived annotations grouped by input molecule

| Molecule | Annotation count | Annotations |
| --- | --- | --- |
| User1 | 0 |  |
