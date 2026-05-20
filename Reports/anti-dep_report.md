# CSANNO readable report

## Metadata

| Field | Value |
| --- | --- |
| Query date/time | 2026-05-20T16:53:54+01:00 |
| Input file | data\anti-dep.sar |
| Molecule count | 5 |
| CSANNO version | CSANNO - (C) 2019/2026 - Andre O. Falcao BioISI - DI/FCUL version 0.5.20260515 |
| ChEMBL server | https://www.ebi.ac.uk |
| ChEMBL request warnings | none |

## CSANNO options

**Command used:**

```bash
csanno.py -in 'data\anti-dep.sar' -sim 0.7 -report AD -digest -glossary
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
| ADRA1A | 5 | 1.0000 |
| DRD1 | 5 | 1.0000 |
| HTR2A | 5 | 1.0000 |
| SLC6A4 | 5 | 1.0000 |
| HRH1 | 4 | 0.8000 |
| HTR2B | 4 | 0.8000 |
| HTR2C | 4 | 0.8000 |
| KCNH2 | 4 | 0.8000 |
| REP | 4 | 0.8000 |
| SIGMAR1 | 4 | 0.8000 |
| SLC6A2 | 4 | 0.8000 |
| SLC6A3 | 4 | 0.8000 |
| ADRA1B | 3 | 0.6000 |
| ADRA2A | 3 | 0.6000 |
| ADRA2B | 3 | 0.6000 |
| TXNRD1 | 3 | 0.6000 |
| ADRA1D | 2 | 0.4000 |
| ADRA2C | 2 | 0.4000 |
| CACNA1C | 2 | 0.4000 |
| CHRM1 | 2 | 0.4000 |
| CHRM2 | 2 | 0.4000 |
| CHRM5 | 2 | 0.4000 |
| CNR1 | 2 | 0.4000 |
| CYP2C19 | 2 | 0.4000 |
| CYP2D6 | 2 | 0.4000 |
| DRD2 | 2 | 0.4000 |
| DRD3 | 2 | 0.4000 |
| EHMT2 | 2 | 0.4000 |
| HRH2 | 2 | 0.4000 |
| HRH3 | 2 | 0.4000 |
| LMNA | 2 | 0.4000 |
| MC3R | 2 | 0.4000 |
| OPRK1 | 2 | 0.4000 |
| PFK | 2 | 0.4000 |
| SLC18A2 | 2 | 0.4000 |
| TMEM97 | 2 | 0.4000 |
| ABCB1 | 1 | 0.2000 |
| ACHE | 1 | 0.2000 |
| ADRB1 | 1 | 0.2000 |
| ALB | 1 | 0.2000 |
| ARSA | 1 | 0.2000 |
| CBX1 | 1 | 0.2000 |
| CCKBR | 1 | 0.2000 |
| CHRM3 | 1 | 0.2000 |
| CHRM4 | 1 | 0.2000 |
| CYP1A2 | 1 | 0.2000 |
| CYP3A4 | 1 | 0.2000 |
| DRD4 | 1 | 0.2000 |
| EDNRA | 1 | 0.2000 |
| FAAH | 1 | 0.2000 |
| GLP1R | 1 | 0.2000 |
| GMNN | 1 | 0.2000 |
| HDAC6 | 1 | 0.2000 |
| HTR1A | 1 | 0.2000 |
| HTR1B | 1 | 0.2000 |
| HTR3A | 1 | 0.2000 |
| HTR6 | 1 | 0.2000 |
| KCNK2 | 1 | 0.2000 |
| KCNK9 | 1 | 0.2000 |
| KMT2A | 1 | 0.2000 |
| MC5R | 1 | 0.2000 |
| OPRM1 | 1 | 0.2000 |
| SLC22A1 | 1 | 0.2000 |
| SLCO1B1 | 1 | 0.2000 |
| SLCO1B3 | 1 | 0.2000 |
| SLCO2B1 | 1 | 0.2000 |
| SMN2 | 1 | 0.2000 |
| TACR1 | 1 | 0.2000 |

### Tier 0 detail - exact ChEMBL matches by input molecule

| Molecule | Annotation count | Annotations |
| --- | --- | --- |
| Sertraline | 35 | ADRA1A, ADRA1B, ADRA2A, ADRA2B, ADRA2C, ADRB1, CACNA1C, CBX1, CHRM1, CHRM2, CHRM4, CHRM5, CYP2C19, CYP2D6, CYP3A4, DRD1, GLP1R, GMNN, HDAC6, HTR2A, HTR2B, HTR2C, KCNH2, KMT2A, MC3R, MC5R, OPRK1, OPRM1, REP, SIGMAR1, SLC18A2, SLC6A2, SLC6A3, SLC6A4, TMEM97 |
| Citalopran | 22 | ADRA1A, ADRA1B, ADRA1D, CHRM2, DRD1, DRD3, EHMT2, HRH1, HRH2, HTR2A, HTR2B, HTR2C, KCNH2, REP, SIGMAR1, SLC22A1, SLC6A2, SLC6A3, SLC6A4, SMN2, TMEM97, TXNRD1 |
| Fluoxetine | 42 | ABCB1, ACHE, ADRA1A, ADRA2A, ADRA2B, CACNA1C, CCKBR, CHRM1, CHRM3, CHRM5, CNR1, CYP1A2, CYP2C19, CYP2D6, DRD1, DRD2, EDNRA, EHMT2, HRH1, HRH2, HRH3, HTR2A, HTR2B, HTR2C, HTR3A, HTR6, KCNH2, KCNK2, KCNK9, LMNA, MC3R, OPRK1, REP, SIGMAR1, SLC18A2, SLC6A2, SLC6A3, SLC6A4, SLCO1B1, SLCO1B3, SLCO2B1, TXNRD1 |
| Tradozone | 23 | ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C, ALB, DRD1, DRD2, DRD3, DRD4, FAAH, HRH1, HTR1A, HTR1B, HTR2A, HTR2B, HTR2C, KCNH2, PFK, SIGMAR1, SLC6A4, TACR1 |
| Venlafaxine | 14 | ADRA1A, ARSA, CNR1, DRD1, HRH1, HRH3, HTR2A, LMNA, PFK, REP, SLC6A2, SLC6A3, SLC6A4, TXNRD1 |

### Tier 1 aggregate — similar ChEMBL analogue pool

| Target | Molecule count | Fraction |
| --- | --- | --- |
| SLC6A4 | 65 | 0.6311 |
| SLC6A2 | 33 | 0.3204 |
| SLC6A3 | 25 | 0.2427 |
| HTR2A | 10 | 0.0971 |
| ADRA1A | 9 | 0.0874 |
| REP | 8 | 0.0777 |
| GMNN | 7 | 0.0680 |
| DRD2 | 6 | 0.0583 |
| HDAC6 | 5 | 0.0485 |
| HRH1 | 5 | 0.0485 |
| HTR6 | 5 | 0.0485 |
| CHRM1 | 4 | 0.0388 |
| GLP1R | 4 | 0.0388 |
| HTR1A | 4 | 0.0388 |
| HTR7 | 4 | 0.0388 |
| KCNH2 | 4 | 0.0388 |
| SIGMAR1 | 4 | 0.0388 |
| TDP1 | 4 | 0.0388 |
| CYP1A2 | 3 | 0.0291 |
| CYP2C19 | 3 | 0.0291 |
| EHMT2 | 3 | 0.0291 |
| LMNA | 3 | 0.0291 |
| NFKB1 | 3 | 0.0291 |
| SLCO1B1 | 3 | 0.0291 |
| SLCO1B3 | 3 | 0.0291 |
| AMPC | 2 | 0.0194 |
| CBX1 | 2 | 0.0194 |
| CYP2D6 | 2 | 0.0194 |
| CYP3A4 | 2 | 0.0194 |
| HTR2C | 2 | 0.0194 |
| MAPK1 | 2 | 0.0194 |
| TXNRD1 | 2 | 0.0194 |
| ALOX15 | 1 | 0.0097 |
| ALOX15B | 1 | 0.0097 |
| APEX1 | 1 | 0.0097 |
| APLNR | 1 | 0.0097 |
| ARSA | 1 | 0.0097 |
| ATXN2 | 1 | 0.0097 |
| BLM | 1 | 0.0097 |
| CX3CR1 | 1 | 0.0097 |
| CYP2C9 | 1 | 0.0097 |
| DRD1 | 1 | 0.0097 |
| FFP | 1 | 0.0097 |
| GBA1 | 1 | 0.0097 |
| GPR183 | 1 | 0.0097 |
| GPR35 | 1 | 0.0097 |
| GPR65 | 1 | 0.0097 |
| HLA-A | 1 | 0.0097 |
| HPGD | 1 | 0.0097 |
| IDH1 | 1 | 0.0097 |
| KAT2A | 1 | 0.0097 |
| MBNL1 | 1 | 0.0097 |
| NFE2L2 | 1 | 0.0097 |
| NS1 | 1 | 0.0097 |
| OPRM1 | 1 | 0.0097 |
| PFK | 1 | 0.0097 |
| POLK | 1 | 0.0097 |
| RGS4 | 1 | 0.0097 |
| RORC | 1 | 0.0097 |
| SMN2 | 1 | 0.0097 |
| TP53 | 1 | 0.0097 |
| TSHR | 1 | 0.0097 |
| USP1 | 1 | 0.0097 |

### Tier 1 detail — analogue-derived annotations grouped by input molecule

| Molecule | Annotation count | Annotations |
| --- | --- | --- |
| Sertraline | 22 | CHRM1, CYP1A2, CYP2C19, CYP2D6, DRD1, EHMT2, FFP, GLP1R, GMNN, HDAC6, IDH1, MAPK1, NFKB1, NS1, REP, SLC6A2, SLC6A3, SLC6A4, SLCO1B1, SLCO1B3, TDP1, TP53 |
| Citalopran | 26 | ADRA1A, ALOX15, CBX1, CHRM1, CYP2C19, CYP2C9, CYP3A4, GMNN, HDAC6, HRH1, HTR2C, HTR6, KCNH2, LMNA, MAPK1, NFKB1, OPRM1, POLK, REP, SLC6A2, SLC6A3, SLC6A4, SLCO1B1, SLCO1B3, TDP1, USP1 |
| Fluoxetine | 40 | ALOX15B, AMPC, APEX1, APLNR, ARSA, ATXN2, BLM, CBX1, CHRM1, CX3CR1, CYP1A2, CYP2C19, CYP2D6, CYP3A4, EHMT2, GBA1, GLP1R, GMNN, GPR183, GPR35, GPR65, HDAC6, HPGD, HTR2C, KAT2A, KCNH2, LMNA, MBNL1, NFKB1, PFK, REP, RGS4, RORC, SLC6A2, SLC6A3, SLC6A4, SMN2, TDP1, TSHR, TXNRD1 |
| Tradozone | 16 | DRD2, EHMT2, GMNN, HDAC6, HLA-A, HRH1, HTR1A, HTR2A, HTR6, HTR7, KCNH2, REP, SIGMAR1, SLC6A2, SLCO1B1, SLCO1B3 |
| Venlafaxine | 8 | ADRA1A, GMNN, HDAC6, NFE2L2, REP, SLC6A2, SLC6A3, SLC6A4 |

## Comprehensive digest

**Tier 0 (direct ChEMBL active evidence).** 5/5 input molecules had at least one active target annotation. The median was 23.0 annotations per molecule, with a maximum of 42.

**Tier 0 promiscuity:** Fluoxetine (42); Sertraline (35); Tradozone (23); Citalopran (22); Venlafaxine (14).

**Tier 0 most common targets:** **ADRA1A** (5 input molecules, 100.0%); **DRD1** (5 input molecules, 100.0%); **HTR2A** (5 input molecules, 100.0%); **SLC6A4** (5 input molecules, 100.0%); **HRH1** (4 input molecules, 80.0%); **HTR2B** (4 input molecules, 80.0%); **HTR2C** (4 input molecules, 80.0%); **KCNH2** (4 input molecules, 80.0%); **REP** (4 input molecules, 80.0%); **SIGMAR1** (4 input molecules, 80.0%).

**Tier 1 (analogue-derived probable target evidence).** 103 similar ChEMBL molecules were evaluated; 85 of them had at least one active target annotation.

Grouped back to the input file, 5/5 input molecules had at least one Tier 1 active target annotation; the median was 22.0 annotations per molecule.

**Tier 1 most common probable targets by input molecule:** **GMNN** (5 input molecules, 100.0%); **HDAC6** (5 input molecules, 100.0%); **REP** (5 input molecules, 100.0%); **SLC6A2** (5 input molecules, 100.0%); **SLC6A3** (4 input molecules, 80.0%); **SLC6A4** (4 input molecules, 80.0%); **CHRM1** (3 input molecules, 60.0%); **CYP2C19** (3 input molecules, 60.0%); **EHMT2** (3 input molecules, 60.0%); **KCNH2** (3 input molecules, 60.0%).

**Tier 1 strongest signals in the analogue pool:** **SLC6A4** (65 similar molecules, 63.1%); **SLC6A2** (33 similar molecules, 32.0%); **SLC6A3** (25 similar molecules, 24.3%); **HTR2A** (10 similar molecules, 9.7%); **ADRA1A** (9 similar molecules, 8.7%); **REP** (8 similar molecules, 7.8%); **GMNN** (7 similar molecules, 6.8%); **DRD2** (6 similar molecules, 5.8%); **HDAC6** (5 similar molecules, 4.9%); **HRH1** (5 similar molecules, 4.9%).

**Tier 1 annotations not seen in Tier 0:** ALOX15, ALOX15B, AMPC, APEX1, APLNR, ATXN2, BLM, CX3CR1, CYP2C9, FFP, GBA1, GPR183, GPR35, GPR65, HLA-A, HPGD, HTR7, IDH1, KAT2A, MAPK1, MBNL1, NFE2L2, NFKB1, NS1, POLK ....

Tier 1 findings are analogue-derived suggestions, not direct activity claims for the input molecules.

## Gene glossary

- [ABCB1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ABCB1) - ATP-dependent translocase ABCB1 (from [UniProt entry P08183](https://www.uniprot.org/uniprotkb/P08183))
- [ACHE](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ACHE) - Acetylcholinesterase (from [UniProt entry P22303](https://www.uniprot.org/uniprotkb/P22303))
- [ADRA1A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ADRA1A) - Alpha-1D adrenergic receptor (from [UniProt entry P25100](https://www.uniprot.org/uniprotkb/P25100))
- [ADRA1B](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ADRA1B) - Alpha-1B adrenergic receptor (from [UniProt entry P35368](https://www.uniprot.org/uniprotkb/P35368))
- [ADRA1D](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ADRA1D) - Alpha-1D adrenergic receptor (from [UniProt entry P25100](https://www.uniprot.org/uniprotkb/P25100))
- [ADRA2A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ADRA2A) - Alpha-2A adrenergic receptor (from [UniProt entry P08913](https://www.uniprot.org/uniprotkb/P08913))
- [ADRA2B](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ADRA2B) - Alpha-2B adrenergic receptor (from [UniProt entry P18089](https://www.uniprot.org/uniprotkb/P18089))
- [ADRA2C](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ADRA2C) - Alpha-2C adrenergic receptor (from [UniProt entry P18825](https://www.uniprot.org/uniprotkb/P18825))
- [ADRB1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ADRB1) - Beta-1 adrenergic receptor (from [UniProt entry P08588](https://www.uniprot.org/uniprotkb/P08588))
- [ALB](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ALB) - Fas-binding factor 1 (from [UniProt entry Q8TES7](https://www.uniprot.org/uniprotkb/Q8TES7))
- [ALOX15](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ALOX15) - Polyunsaturated fatty acid lipoxygenase ALOX15 (from [UniProt entry P16050](https://www.uniprot.org/uniprotkb/P16050))
- [ALOX15B](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ALOX15B) - Polyunsaturated fatty acid lipoxygenase ALOX15B (from [UniProt entry O15296](https://www.uniprot.org/uniprotkb/O15296))
- [AMPC](https://www.genecards.org/cgi-bin/carddisp.pl?gene=AMPC) - Beta-lactamase (from [UniProt entry P00811](https://www.uniprot.org/uniprotkb/P00811))
- [APEX1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=APEX1) - DNA repair nuclease/redox regulator APEX1 (from [UniProt entry P27695](https://www.uniprot.org/uniprotkb/P27695))
- [APLNR](https://www.genecards.org/cgi-bin/carddisp.pl?gene=APLNR) - Apelin receptor (from [UniProt entry P35414](https://www.uniprot.org/uniprotkb/P35414))
- [ARSA](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ARSA) - ATPase GET3 (from [UniProt entry O43681](https://www.uniprot.org/uniprotkb/O43681))
- [ATXN2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ATXN2) - Ataxin-2 (from [UniProt entry Q99700](https://www.uniprot.org/uniprotkb/Q99700))
- [BLM](https://www.genecards.org/cgi-bin/carddisp.pl?gene=BLM) - RecQ-like DNA helicase BLM (from [UniProt entry P54132](https://www.uniprot.org/uniprotkb/P54132))
- [CACNA1C](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CACNA1C) - Voltage-dependent L-type calcium channel subunit alpha-1C (from [UniProt entry Q13936](https://www.uniprot.org/uniprotkb/Q13936))
- [CBX1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CBX1) - Chromobox protein homolog 1 (from [UniProt entry P83916](https://www.uniprot.org/uniprotkb/P83916))
- [CCKBR](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CCKBR) - Gastrin/cholecystokinin type B receptor (from [UniProt entry P32239](https://www.uniprot.org/uniprotkb/P32239))
- [CHRM1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CHRM1) - Muscarinic acetylcholine receptor M1 (from [UniProt entry P11229](https://www.uniprot.org/uniprotkb/P11229))
- [CHRM2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CHRM2) - Muscarinic acetylcholine receptor M2 (from [UniProt entry P08172](https://www.uniprot.org/uniprotkb/P08172))
- [CHRM3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CHRM3) - Muscarinic acetylcholine receptor M3 (from [UniProt entry P20309](https://www.uniprot.org/uniprotkb/P20309))
- [CHRM4](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CHRM4) - Muscarinic acetylcholine receptor M4 (from [UniProt entry P08173](https://www.uniprot.org/uniprotkb/P08173))
- [CHRM5](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CHRM5) - Muscarinic acetylcholine receptor M5 (from [UniProt entry P08912](https://www.uniprot.org/uniprotkb/P08912))
- [CNR1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CNR1) - Cannabinoid receptor 1 (from [UniProt entry P21554](https://www.uniprot.org/uniprotkb/P21554))
- [CX3CR1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CX3CR1) - CX3C chemokine receptor 1 (from [UniProt entry P49238](https://www.uniprot.org/uniprotkb/P49238))
- [CYP1A2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYP1A2) - Cytochrome P450 1A2 (from [UniProt entry P05177](https://www.uniprot.org/uniprotkb/P05177))
- [CYP2C19](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYP2C19) - Cytochrome P450 2C19 (from [UniProt entry P33261](https://www.uniprot.org/uniprotkb/P33261))
- [CYP2C9](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYP2C9) - Cytochrome P450 2C9 (from [UniProt entry P11712](https://www.uniprot.org/uniprotkb/P11712))
- [CYP2D6](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYP2D6) - Cytochrome P450 2D6 (from [UniProt entry P10635](https://www.uniprot.org/uniprotkb/P10635))
- [CYP3A4](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYP3A4) - Cytochrome P450 3A4 (from [UniProt entry P08684](https://www.uniprot.org/uniprotkb/P08684))
- [DRD1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=DRD1) - D(1A) dopamine receptor (from [UniProt entry P21728](https://www.uniprot.org/uniprotkb/P21728))
- [DRD2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=DRD2) - D(2) dopamine receptor (from [UniProt entry P14416](https://www.uniprot.org/uniprotkb/P14416))
- [DRD3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=DRD3) - D(3) dopamine receptor (from [UniProt entry P35462](https://www.uniprot.org/uniprotkb/P35462))
- [DRD4](https://www.genecards.org/cgi-bin/carddisp.pl?gene=DRD4) - D(4) dopamine receptor (from [UniProt entry P21917](https://www.uniprot.org/uniprotkb/P21917))
- [EDNRA](https://www.genecards.org/cgi-bin/carddisp.pl?gene=EDNRA) - Endothelin-1 receptor (from [UniProt entry P25101](https://www.uniprot.org/uniprotkb/P25101))
- [EHMT2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=EHMT2) - Histone-lysine N-methyltransferase EHMT2 (from [UniProt entry Q96KQ7](https://www.uniprot.org/uniprotkb/Q96KQ7))
- [FAAH](https://www.genecards.org/cgi-bin/carddisp.pl?gene=FAAH) - Fatty acid 2-hydroxylase (from [UniProt entry Q7L5A8](https://www.uniprot.org/uniprotkb/Q7L5A8))
- [FFP](https://www.genecards.org/cgi-bin/carddisp.pl?gene=FFP) - 4'-phosphopantetheinyl transferase ffp (from [UniProt entry Q9F4F7](https://www.uniprot.org/uniprotkb/Q9F4F7))
- [GBA1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GBA1) - Lysosomal acid glucosylceramidase (from [UniProt entry P04062](https://www.uniprot.org/uniprotkb/P04062))
- [GLP1R](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GLP1R) - Glucagon-like peptide 1 receptor (from [UniProt entry P43220](https://www.uniprot.org/uniprotkb/P43220))
- [GMNN](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GMNN) - Geminin (from [UniProt entry O75496](https://www.uniprot.org/uniprotkb/O75496))
- [GPR183](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GPR183) - G-protein coupled receptor 183 (from [UniProt entry P32249](https://www.uniprot.org/uniprotkb/P32249))
- [GPR35](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GPR35) - G-protein coupled receptor 35 (from [UniProt entry Q9HC97](https://www.uniprot.org/uniprotkb/Q9HC97))
- [GPR65](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GPR65) - G-protein coupled receptor 65 (from [UniProt entry Q8IYL9](https://www.uniprot.org/uniprotkb/Q8IYL9))
- [HDAC6](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HDAC6) - Protein deacetylase HDAC6 (from [UniProt entry Q9UBN7](https://www.uniprot.org/uniprotkb/Q9UBN7))
- [HLA-A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HLA-A) - HLA class I histocompatibility antigen, A alpha chain (from [UniProt entry P04439](https://www.uniprot.org/uniprotkb/P04439))
- [HPGD](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HPGD) - 15-hydroxyprostaglandin dehydrogenase [NAD(+)] (from [UniProt entry P15428](https://www.uniprot.org/uniprotkb/P15428))
- [HRH1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HRH1) - Histamine H1 receptor (from [UniProt entry P35367](https://www.uniprot.org/uniprotkb/P35367))
- [HRH2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HRH2) - Histamine H2 receptor (from [UniProt entry P25021](https://www.uniprot.org/uniprotkb/P25021))
- [HRH3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HRH3) - Histamine H3 receptor (from [UniProt entry Q9Y5N1](https://www.uniprot.org/uniprotkb/Q9Y5N1))
- [HTR1A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR1A) - 5-hydroxytryptamine receptor 1A (from [UniProt entry P08908](https://www.uniprot.org/uniprotkb/P08908))
- [HTR1B](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR1B) - 5-hydroxytryptamine receptor 1B (from [UniProt entry P28222](https://www.uniprot.org/uniprotkb/P28222))
- [HTR2A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR2A) - 5-hydroxytryptamine receptor 2A (from [UniProt entry P28223](https://www.uniprot.org/uniprotkb/P28223))
- [HTR2B](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR2B) - 5-hydroxytryptamine receptor 2B (from [UniProt entry P41595](https://www.uniprot.org/uniprotkb/P41595))
- [HTR2C](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR2C) - 5-hydroxytryptamine receptor 2C (from [UniProt entry P28335](https://www.uniprot.org/uniprotkb/P28335))
- [HTR3A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR3A) - 5-hydroxytryptamine receptor 3A (from [UniProt entry P46098](https://www.uniprot.org/uniprotkb/P46098))
- [HTR6](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR6) - 5-hydroxytryptamine receptor 6 (from [UniProt entry P50406](https://www.uniprot.org/uniprotkb/P50406))
- [HTR7](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR7) - 5-hydroxytryptamine receptor 7 (from [UniProt entry P34969](https://www.uniprot.org/uniprotkb/P34969))
- [IDH1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=IDH1) - Isocitrate dehydrogenase [NADP] cytoplasmic (from [UniProt entry O75874](https://www.uniprot.org/uniprotkb/O75874))
- [KAT2A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=KAT2A) - Histone acetyltransferase KAT2A (from [UniProt entry Q92830](https://www.uniprot.org/uniprotkb/Q92830))
- [KCNH2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=KCNH2) - Voltage-gated inwardly rectifying potassium channel KCNH2 (from [UniProt entry Q12809](https://www.uniprot.org/uniprotkb/Q12809))
- [KCNK2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=KCNK2) - Potassium channel subfamily K member 2 (from [UniProt entry O95069](https://www.uniprot.org/uniprotkb/O95069))
- [KCNK9](https://www.genecards.org/cgi-bin/carddisp.pl?gene=KCNK9) - Potassium channel subfamily K member 9 (from [UniProt entry Q9NPC2](https://www.uniprot.org/uniprotkb/Q9NPC2))
- [KMT2A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=KMT2A) - Histone-lysine N-methyltransferase 2A (from [UniProt entry Q03164](https://www.uniprot.org/uniprotkb/Q03164))
- [LMNA](https://www.genecards.org/cgi-bin/carddisp.pl?gene=LMNA) - Prelamin-A/C (from [UniProt entry P02545](https://www.uniprot.org/uniprotkb/P02545))
- [MAPK1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MAPK1) - Mitogen-activated protein kinase 1 (from [UniProt entry P28482](https://www.uniprot.org/uniprotkb/P28482))
- [MBNL1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MBNL1) - Muscleblind-like protein 1 (from [UniProt entry Q9NR56](https://www.uniprot.org/uniprotkb/Q9NR56))
- [MC3R](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MC3R) - Melanocortin receptor 3 (from [UniProt entry P41968](https://www.uniprot.org/uniprotkb/P41968))
- [MC5R](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MC5R) - Melanocortin receptor 5 (from [UniProt entry P33032](https://www.uniprot.org/uniprotkb/P33032))
- [NFE2L2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NFE2L2) - Nuclear factor erythroid 2-related factor 2 (from [UniProt entry Q16236](https://www.uniprot.org/uniprotkb/Q16236))
- [NFKB1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NFKB1) - Nuclear factor NF-kappa-B p105 subunit (from [UniProt entry P19838](https://www.uniprot.org/uniprotkb/P19838))
- [NS1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NS1) - Influenza virus NS1A-binding protein (from [UniProt entry Q9Y6Y0](https://www.uniprot.org/uniprotkb/Q9Y6Y0))
- [OPRK1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=OPRK1) - Kappa-type opioid receptor (from [UniProt entry P41145](https://www.uniprot.org/uniprotkb/P41145))
- [OPRM1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=OPRM1) - Mu-type opioid receptor (from [UniProt entry P35372](https://www.uniprot.org/uniprotkb/P35372))
- [PFK](https://www.genecards.org/cgi-bin/carddisp.pl?gene=PFK) - ATP-dependent 6-phosphofructokinase (from [UniProt entry P00512](https://www.uniprot.org/uniprotkb/P00512))
- [POLK](https://www.genecards.org/cgi-bin/carddisp.pl?gene=POLK) - DNA polymerase kappa (from [UniProt entry Q9UBT6](https://www.uniprot.org/uniprotkb/Q9UBT6))
- [REP](https://www.genecards.org/cgi-bin/carddisp.pl?gene=REP) - Replicase polyprotein 1ab (from [UniProt entry P0C6W7](https://www.uniprot.org/uniprotkb/P0C6W7))
- [RGS4](https://www.genecards.org/cgi-bin/carddisp.pl?gene=RGS4) - Regulator of G-protein signaling 4 (from [UniProt entry P49798](https://www.uniprot.org/uniprotkb/P49798))
- [RORC](https://www.genecards.org/cgi-bin/carddisp.pl?gene=RORC) - Nuclear receptor ROR-gamma (from [UniProt entry P51449](https://www.uniprot.org/uniprotkb/P51449))
- [SIGMAR1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SIGMAR1) - Sigma non-opioid intracellular receptor 1 (from [UniProt entry Q99720](https://www.uniprot.org/uniprotkb/Q99720))
- [SLC18A2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC18A2) - Synaptic vesicular amine transporter (from [UniProt entry Q05940](https://www.uniprot.org/uniprotkb/Q05940))
- [SLC22A1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC22A1) - Solute carrier family 22 member 1 (from [UniProt entry O15245](https://www.uniprot.org/uniprotkb/O15245))
- [SLC6A2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC6A2) - Sodium-dependent noradrenaline transporter (from [UniProt entry P23975](https://www.uniprot.org/uniprotkb/P23975))
- [SLC6A3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC6A3) - Sodium-dependent dopamine transporter (from [UniProt entry Q01959](https://www.uniprot.org/uniprotkb/Q01959))
- [SLC6A4](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC6A4) - Sodium-dependent serotonin transporter (from [UniProt entry P31645](https://www.uniprot.org/uniprotkb/P31645))
- [SLCO1B1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLCO1B1) - Solute carrier organic anion transporter family member 1B1 (from [UniProt entry Q9Y6L6](https://www.uniprot.org/uniprotkb/Q9Y6L6))
- [SLCO1B3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLCO1B3) - Solute carrier organic anion transporter family member 1B3 (from [UniProt entry Q9NPD5](https://www.uniprot.org/uniprotkb/Q9NPD5))
- [SLCO2B1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLCO2B1) - Solute carrier organic anion transporter family member 2B1 (from [UniProt entry O94956](https://www.uniprot.org/uniprotkb/O94956))
- [SMN2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SMN2) - Survival motor neuron protein (from [UniProt entry Q16637](https://www.uniprot.org/uniprotkb/Q16637))
- [TACR1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TACR1) - Substance-P receptor (from [UniProt entry P25103](https://www.uniprot.org/uniprotkb/P25103))
- [TDP1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TDP1) - Tyrosyl-DNA phosphodiesterase 1 (from [UniProt entry Q9NUW8](https://www.uniprot.org/uniprotkb/Q9NUW8))
- [TMEM97](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TMEM97) - Sigma intracellular receptor 2 (from [UniProt entry Q5BJF2](https://www.uniprot.org/uniprotkb/Q5BJF2))
- [TP53](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TP53) - Cellular tumor antigen p53 (from [UniProt entry P04637](https://www.uniprot.org/uniprotkb/P04637))
- [TSHR](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TSHR) - Thyrotropin receptor (from [UniProt entry P16473](https://www.uniprot.org/uniprotkb/P16473))
- [TXNRD1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TXNRD1) - Thioredoxin reductase 1, cytoplasmic (from [UniProt entry Q16881](https://www.uniprot.org/uniprotkb/Q16881))
- [USP1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=USP1) - Ubiquitin carboxyl-terminal hydrolase 1 (from [UniProt entry O94782](https://www.uniprot.org/uniprotkb/O94782))
