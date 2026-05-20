# CSANNO readable report

## Metadata

| Field | Value |
| --- | --- |
| Query date/time | 2026-05-20T11:29:17+01:00 |
| Input file | data\anticonvulsants.sar |
| Molecule count | 12 |
| CSANNO version | CSANNO - (C) 2019/2026 - Andre O. Falcao BioISI - DI/FCUL version 0.5.20260515 |
| ChEMBL server | https://www.ebi.ac.uk |
| ChEMBL request warnings | none |

## CSANNO options

**Command used:**

```bash
csanno.py -in 'data\anticonvulsants.sar' -out 'data\anticonv-A' -sim 0.7 -digest -glossary -report AD
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
| REP | 10 | 0.8333 |
| SLCO1B3 | 8 | 0.6667 |
| SLCO1B1 | 7 | 0.5833 |
| SAMHD1 | 6 | 0.5000 |
| EHMT2 | 5 | 0.4167 |
| GMNN | 5 | 0.4167 |
| SLCO2B1 | 5 | 0.4167 |
| APEX1 | 3 | 0.2500 |
| CA1 | 3 | 0.2500 |
| CA12 | 3 | 0.2500 |
| CA13 | 3 | 0.2500 |
| CA14 | 3 | 0.2500 |
| CA15 | 3 | 0.2500 |
| CA2 | 3 | 0.2500 |
| CA4 | 3 | 0.2500 |
| CA5A | 3 | 0.2500 |
| CA5B | 3 | 0.2500 |
| CA6 | 3 | 0.2500 |
| CA7 | 3 | 0.2500 |
| CA9 | 3 | 0.2500 |
| CHRM1 | 3 | 0.2500 |
| HDAC6 | 3 | 0.2500 |
| LMNA | 3 | 0.2500 |
| RGS4 | 3 | 0.2500 |
| SLC22A1 | 3 | 0.2500 |
| THRB | 3 | 0.2500 |
| TXNRD1 | 3 | 0.2500 |
| ADRA1A | 2 | 0.1667 |
| CA | 2 | 0.1667 |
| CA3 | 2 | 0.1667 |
| CACNA2D1 | 2 | 0.1667 |
| CAN2 | 2 | 0.1667 |
| CYNT | 2 | 0.1667 |
| GALC | 2 | 0.1667 |
| HRH1 | 2 | 0.1667 |
| HTR2A | 2 | 0.1667 |
| MDH1 | 2 | 0.1667 |
| MTCA2 | 2 | 0.1667 |
| NCE103 | 2 | 0.1667 |
| P2RX4 | 2 | 0.1667 |
| SCN2A | 2 | 0.1667 |
| SCN3A | 2 | 0.1667 |
| SCN9A | 2 | 0.1667 |
| 1272966 | 1 | 0.0833 |
| ALDH1A1 | 1 | 0.0833 |
| AR | 1 | 0.0833 |
| ARSA | 1 | 0.0833 |
| ATAD5 | 1 | 0.0833 |
| BLM | 1 | 0.0833 |
| CACNA2D2 | 1 | 0.0833 |
| CHRM3 | 1 | 0.0833 |
| CYP2C9 | 1 | 0.0833 |
| CYP2D6 | 1 | 0.0833 |
| CYP3A4 | 1 | 0.0833 |
| DRD1 | 1 | 0.0833 |
| DRD3 | 1 | 0.0833 |
| F3 | 1 | 0.0833 |
| FAAH | 1 | 0.0833 |
| FEN1 | 1 | 0.0833 |
| GSK3B | 1 | 0.0833 |
| IDH1 | 1 | 0.0833 |
| IDO2 | 1 | 0.0833 |
| IMPA1 | 1 | 0.0833 |
| KCNH2 | 1 | 0.0833 |
| MAOB | 1 | 0.0833 |
| MAPK1 | 1 | 0.0833 |
| MTCA1 | 1 | 0.0833 |
| NFO | 1 | 0.0833 |
| P2RX7 | 1 | 0.0833 |
| PDE3A | 1 | 0.0833 |
| PFK | 1 | 0.0833 |
| PMP22 | 1 | 0.0833 |
| PTH1R | 1 | 0.0833 |
| RORC | 1 | 0.0833 |
| SCN10A | 1 | 0.0833 |
| SCN1A | 1 | 0.0833 |
| SCN4A | 1 | 0.0833 |
| SCN5A | 1 | 0.0833 |
| SCN8A | 1 | 0.0833 |
| SLC22A6 | 1 | 0.0833 |
| SLC6A3 | 1 | 0.0833 |
| SMN2 | 1 | 0.0833 |
| TDP1 | 1 | 0.0833 |
| TP53 | 1 | 0.0833 |
| TSHR | 1 | 0.0833 |
| UGT1A9 | 1 | 0.0833 |
| UGT2B15 | 1 | 0.0833 |
| VPR | 1 | 0.0833 |

### Tier 0 detail - exact ChEMBL matches by input molecule

| Molecule | Annotation count | Annotations |
| --- | --- | --- |
| Carbamazepine | 18 | AR, EHMT2, GMNN, IDH1, LMNA, MAPK1, MDH1, P2RX4, REP, SAMHD1, SLC22A1, SLCO1B1, SLCO1B3, SLCO2B1, THRB, TP53, TSHR, VPR |
| Valproic acid | 22 | ADRA1A, ATAD5, CHRM1, CHRM3, DRD3, F3, GSK3B, HRH1, HTR2A, RGS4, SAMHD1, SLC22A1, SLC22A6, SLC6A3, SLCO1B1, SLCO1B3, SLCO2B1, TDP1, THRB, TXNRD1, UGT1A9, UGT2B15 |
| Lamotrigine | 23 | APEX1, CHRM1, CYP2D6, DRD1, EHMT2, GALC, GMNN, MDH1, PDE3A, REP, SAMHD1, SCN10A, SCN1A, SCN2A, SCN3A, SCN4A, SCN5A, SCN8A, SCN9A, SLCO1B1, SLCO1B3, SLCO2B1, TXNRD1 |
| Levetiracetam | 4 | EHMT2, IDO2, REP, SAMHD1 |
| Phenytoin | 25 | ARSA, BLM, CHRM1, CYP2C9, EHMT2, FAAH, FEN1, GMNN, HRH1, HTR2A, IMPA1, KCNH2, LMNA, NFO, PTH1R, REP, RGS4, RORC, SAMHD1, SCN2A, SLCO1B1, SLCO1B3, SLCO2B1, THRB, TXNRD1 |
| Topiramate | 23 | APEX1, CA, CA1, CA12, CA13, CA14, CA15, CA2, CA4, CA5A, CA5B, CA6, CA7, CA9, CAN2, CYNT, HDAC6, MTCA1, MTCA2, NCE103, REP, SLCO1B1, SLCO1B3 |
| Gabapentin | 8 | ALDH1A1, CACNA2D1, EHMT2, PFK, REP, RGS4, SLC22A1, SMN2 |
| Pregabalin | 9 | ADRA1A, CACNA2D1, CACNA2D2, CYP3A4, HDAC6, P2RX7, REP, SLCO1B1, SLCO1B3 |
| Phenobarbital | 3 | GMNN, SLCO1B3, SLCO2B1 |
| Oxcarbazepine | 7 | GALC, HDAC6, P2RX4, PMP22, REP, SLCO1B1, SLCO1B3 |
| Lacosamide | 16 | CA1, CA12, CA13, CA14, CA15, CA2, CA3, CA4, CA5A, CA5B, CA6, CA7, CA9, REP, SCN3A, SCN9A |
| Zonisamide | 25 | 1272966, APEX1, CA, CA1, CA12, CA13, CA14, CA15, CA2, CA3, CA4, CA5A, CA5B, CA6, CA7, CA9, CAN2, CYNT, GMNN, LMNA, MAOB, MTCA2, NCE103, REP, SAMHD1 |

### Tier 1 aggregate — similar ChEMBL analogue pool

| Target | Molecule count | Fraction |
| --- | --- | --- |
| CA2 | 5 | 0.0368 |
| P2RX4 | 4 | 0.0294 |
| REP | 3 | 0.0221 |
| BLM | 2 | 0.0147 |
| SLCO1B1 | 2 | 0.0147 |
| SLCO1B3 | 2 | 0.0147 |
| AKR1B1 | 1 | 0.0074 |
| APOBEC3F | 1 | 0.0074 |
| ATXN2 | 1 | 0.0074 |
| CACNA2D1 | 1 | 0.0074 |
| CHRM1 | 1 | 0.0074 |
| CYP3A4 | 1 | 0.0074 |
| EHMT2 | 1 | 0.0074 |
| GLA | 1 | 0.0074 |
| GMNN | 1 | 0.0074 |
| GNAS | 1 | 0.0074 |
| KAT2A | 1 | 0.0074 |
| KMT2A | 1 | 0.0074 |
| LMNA | 1 | 0.0074 |
| LTA4H | 1 | 0.0074 |
| MMP2 | 1 | 0.0074 |
| MMP3 | 1 | 0.0074 |
| MMP9 | 1 | 0.0074 |
| NFE2L2 | 1 | 0.0074 |
| NFKB1 | 1 | 0.0074 |
| NFO | 1 | 0.0074 |
| NPSR1 | 1 | 0.0074 |
| NS1 | 1 | 0.0074 |
| SMN2 | 1 | 0.0074 |
| THRB | 1 | 0.0074 |

### Tier 1 detail — analogue-derived annotations grouped by input molecule

| Molecule | Annotation count | Annotations |
| --- | --- | --- |
| Carbamazepine | 1 | P2RX4 |
| Valproic acid | 10 | BLM, CHRM1, CYP3A4, KAT2A, LTA4H, NFKB1, NPSR1, REP, SLCO1B1, SLCO1B3 |
| Lamotrigine | 0 |  |
| Levetiracetam | 1 | REP |
| Phenytoin | 11 | AKR1B1, EHMT2, GLA, GMNN, KMT2A, NFE2L2, NS1, REP, SLCO1B1, SLCO1B3, THRB |
| Topiramate | 1 | CA2 |
| Gabapentin | 1 | CACNA2D1 |
| Pregabalin | 1 | LMNA |
| Phenobarbital | 6 | APOBEC3F, ATXN2, GNAS, MMP2, MMP3, MMP9 |
| Oxcarbazepine | 0 |  |
| Lacosamide | 0 |  |
| Zonisamide | 3 | BLM, NFO, SMN2 |

## Comprehensive digest

**Tier 0 (direct ChEMBL active evidence).** 12/12 input molecules had at least one active target annotation. The median was 17.0 annotations per molecule, with a maximum of 25.

**Tier 0 promiscuity:** Phenytoin (25); Zonisamide (25); Lamotrigine (23); Topiramate (23); Valproic acid (22).

**Tier 0 most common targets:** **REP** (10 input molecules, 83.3%); **SLCO1B3** (8 input molecules, 66.7%); **SLCO1B1** (7 input molecules, 58.3%); **SAMHD1** (6 input molecules, 50.0%); **EHMT2** (5 input molecules, 41.7%); **GMNN** (5 input molecules, 41.7%); **SLCO2B1** (5 input molecules, 41.7%); **APEX1** (3 input molecules, 25.0%); **CA1** (3 input molecules, 25.0%); **CA12** (3 input molecules, 25.0%).

**Tier 1 (analogue-derived probable target evidence).** 136 similar ChEMBL molecules were evaluated; 18 of them had at least one active target annotation.

Grouped back to the input file, 9/12 input molecules had at least one Tier 1 active target annotation; the median was 1.0 annotations per molecule.

**Tier 1 most common probable targets by input molecule:** **REP** (3 input molecules, 25.0%); **BLM** (2 input molecules, 16.7%); **SLCO1B1** (2 input molecules, 16.7%); **SLCO1B3** (2 input molecules, 16.7%); **AKR1B1** (1 input molecules, 8.3%); **APOBEC3F** (1 input molecules, 8.3%); **ATXN2** (1 input molecules, 8.3%); **CA2** (1 input molecules, 8.3%); **CACNA2D1** (1 input molecules, 8.3%); **CHRM1** (1 input molecules, 8.3%).

**Tier 1 strongest signals in the analogue pool:** **CA2** (5 similar molecules, 3.7%); **P2RX4** (4 similar molecules, 2.9%); **REP** (3 similar molecules, 2.2%); **BLM** (2 similar molecules, 1.5%); **SLCO1B1** (2 similar molecules, 1.5%); **SLCO1B3** (2 similar molecules, 1.5%); **AKR1B1** (1 similar molecules, 0.7%); **APOBEC3F** (1 similar molecules, 0.7%); **ATXN2** (1 similar molecules, 0.7%); **CACNA2D1** (1 similar molecules, 0.7%).

**Tier 1 annotations not seen in Tier 0:** AKR1B1, APOBEC3F, ATXN2, GLA, GNAS, KAT2A, KMT2A, LTA4H, MMP2, MMP3, MMP9, NFE2L2, NFKB1, NPSR1, NS1.

Tier 1 findings are analogue-derived suggestions, not direct activity claims for the input molecules.

## Gene glossary

- [ADRA1A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ADRA1A) - Alpha-1D adrenergic receptor (from [UniProt entry P25100](https://www.uniprot.org/uniprotkb/P25100))
- [AKR1B1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=AKR1B1) - Aldo-keto reductase family 1 member B1 (from [UniProt entry P15121](https://www.uniprot.org/uniprotkb/P15121))
- [ALDH1A1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ALDH1A1) - Aldehyde dehydrogenase 1A1 (from [UniProt entry P00352](https://www.uniprot.org/uniprotkb/P00352))
- [APEX1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=APEX1) - DNA repair nuclease/redox regulator APEX1 (from [UniProt entry P27695](https://www.uniprot.org/uniprotkb/P27695))
- [APOBEC3F](https://www.genecards.org/cgi-bin/carddisp.pl?gene=APOBEC3F) - DNA dC->dU-editing enzyme APOBEC-3F (from [UniProt entry Q8IUX4](https://www.uniprot.org/uniprotkb/Q8IUX4))
- [AR](https://www.genecards.org/cgi-bin/carddisp.pl?gene=AR) - Androgen receptor (from [UniProt entry P10275](https://www.uniprot.org/uniprotkb/P10275))
- [ARSA](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ARSA) - ATPase GET3 (from [UniProt entry O43681](https://www.uniprot.org/uniprotkb/O43681))
- [ATAD5](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ATAD5) - ATPase family AAA domain-containing protein 5 (from [UniProt entry Q96QE3](https://www.uniprot.org/uniprotkb/Q96QE3))
- [ATXN2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ATXN2) - Ataxin-2 (from [UniProt entry Q99700](https://www.uniprot.org/uniprotkb/Q99700))
- [BLM](https://www.genecards.org/cgi-bin/carddisp.pl?gene=BLM) - RecQ-like DNA helicase BLM (from [UniProt entry P54132](https://www.uniprot.org/uniprotkb/P54132))
- [CA](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA) - Carbonic anhydrase (from [UniProt entry A0A1S6KKC9](https://www.uniprot.org/uniprotkb/A0A1S6KKC9))
- [CA1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA1) - Carbonic anhydrase 1 (from [UniProt entry P00915](https://www.uniprot.org/uniprotkb/P00915))
- [CA12](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA12) - Carbonic anhydrase 12 (from [UniProt entry O43570](https://www.uniprot.org/uniprotkb/O43570))
- [CA13](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA13) - Carbonic anhydrase 13 (from [UniProt entry Q8N1Q1](https://www.uniprot.org/uniprotkb/Q8N1Q1))
- [CA14](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA14) - Carbonic anhydrase 14 (from [UniProt entry Q9ULX7](https://www.uniprot.org/uniprotkb/Q9ULX7))
- [CA15](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA15) - Carbonic anhydrase 15 (from [UniProt entry Q99N23](https://www.uniprot.org/uniprotkb/Q99N23))
- [CA2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA2) - Carbonic anhydrase 2 (from [UniProt entry P00918](https://www.uniprot.org/uniprotkb/P00918))
- [CA3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA3) - Carbonic anhydrase 3 (from [UniProt entry P07451](https://www.uniprot.org/uniprotkb/P07451))
- [CA4](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA4) - Carbonic anhydrase 4 (from [UniProt entry P22748](https://www.uniprot.org/uniprotkb/P22748))
- [CA5A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA5A) - Carbonic anhydrase 5A, mitochondrial (from [UniProt entry P35218](https://www.uniprot.org/uniprotkb/P35218))
- [CA5B](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA5B) - Carbonic anhydrase 5B, mitochondrial (from [UniProt entry Q9Y2D0](https://www.uniprot.org/uniprotkb/Q9Y2D0))
- [CA6](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA6) - Carbonic anhydrase 6 (from [UniProt entry P23280](https://www.uniprot.org/uniprotkb/P23280))
- [CA7](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA7) - Carbonic anhydrase 7 (from [UniProt entry P43166](https://www.uniprot.org/uniprotkb/P43166))
- [CA9](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CA9) - Carbonic anhydrase 9 (from [UniProt entry Q16790](https://www.uniprot.org/uniprotkb/Q16790))
- [CACNA2D1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CACNA2D1) - Voltage-dependent calcium channel subunit alpha-2/delta-1 (from [UniProt entry P54289](https://www.uniprot.org/uniprotkb/P54289))
- [CACNA2D2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CACNA2D2) - Voltage-dependent calcium channel subunit alpha-2/delta-2 (from [UniProt entry Q9NY47](https://www.uniprot.org/uniprotkb/Q9NY47))
- [CAN2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CAN2) - Staphylococcal-like nuclease CAN2 (from [UniProt entry F4IH31](https://www.uniprot.org/uniprotkb/F4IH31))
- [CHRM1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CHRM1) - Muscarinic acetylcholine receptor M1 (from [UniProt entry P11229](https://www.uniprot.org/uniprotkb/P11229))
- [CHRM3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CHRM3) - Muscarinic acetylcholine receptor M3 (from [UniProt entry P20309](https://www.uniprot.org/uniprotkb/P20309))
- [CYNT](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYNT) - Carbonic anhydrase (from [UniProt entry Q9ZN54](https://www.uniprot.org/uniprotkb/Q9ZN54))
- [CYP2C9](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYP2C9) - Cytochrome P450 2C9 (from [UniProt entry P11712](https://www.uniprot.org/uniprotkb/P11712))
- [CYP2D6](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYP2D6) - Cytochrome P450 2D6 (from [UniProt entry P10635](https://www.uniprot.org/uniprotkb/P10635))
- [CYP3A4](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CYP3A4) - Cytochrome P450 3A4 (from [UniProt entry P08684](https://www.uniprot.org/uniprotkb/P08684))
- [DRD1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=DRD1) - D(1A) dopamine receptor (from [UniProt entry P21728](https://www.uniprot.org/uniprotkb/P21728))
- [DRD3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=DRD3) - D(3) dopamine receptor (from [UniProt entry P35462](https://www.uniprot.org/uniprotkb/P35462))
- [EHMT2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=EHMT2) - Histone-lysine N-methyltransferase EHMT2 (from [UniProt entry Q96KQ7](https://www.uniprot.org/uniprotkb/Q96KQ7))
- [F3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=F3) - Tissue factor (from [UniProt entry P13726](https://www.uniprot.org/uniprotkb/P13726))
- [FAAH](https://www.genecards.org/cgi-bin/carddisp.pl?gene=FAAH) - Fatty acid 2-hydroxylase (from [UniProt entry Q7L5A8](https://www.uniprot.org/uniprotkb/Q7L5A8))
- [FEN1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=FEN1) - Flap endonuclease 1 (from [UniProt entry P39748](https://www.uniprot.org/uniprotkb/P39748))
- [GALC](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GALC) - Galactocerebrosidase (from [UniProt entry P54803](https://www.uniprot.org/uniprotkb/P54803))
- [GLA](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GLA) - N-acetyltransferase 8 (from [UniProt entry Q9UHE5](https://www.uniprot.org/uniprotkb/Q9UHE5))
- [GMNN](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GMNN) - Geminin (from [UniProt entry O75496](https://www.uniprot.org/uniprotkb/O75496))
- [GNAS](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GNAS) - Protein ALEX (from [UniProt entry P84996](https://www.uniprot.org/uniprotkb/P84996))
- [GSK3B](https://www.genecards.org/cgi-bin/carddisp.pl?gene=GSK3B) - Glycogen synthase kinase-3 beta (from [UniProt entry P49841](https://www.uniprot.org/uniprotkb/P49841))
- [HDAC6](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HDAC6) - Protein deacetylase HDAC6 (from [UniProt entry Q9UBN7](https://www.uniprot.org/uniprotkb/Q9UBN7))
- [HRH1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HRH1) - Histamine H1 receptor (from [UniProt entry P35367](https://www.uniprot.org/uniprotkb/P35367))
- [HTR2A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR2A) - 5-hydroxytryptamine receptor 2A (from [UniProt entry P28223](https://www.uniprot.org/uniprotkb/P28223))
- [IDH1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=IDH1) - Isocitrate dehydrogenase [NADP] cytoplasmic (from [UniProt entry O75874](https://www.uniprot.org/uniprotkb/O75874))
- [IDO2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=IDO2) - Indoleamine 2,3-dioxygenase 2 (from [UniProt entry Q6ZQW0](https://www.uniprot.org/uniprotkb/Q6ZQW0))
- [IMPA1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=IMPA1) - Inositol monophosphatase 1 (from [UniProt entry P29218](https://www.uniprot.org/uniprotkb/P29218))
- [KAT2A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=KAT2A) - Histone acetyltransferase KAT2A (from [UniProt entry Q92830](https://www.uniprot.org/uniprotkb/Q92830))
- [KCNH2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=KCNH2) - Voltage-gated inwardly rectifying potassium channel KCNH2 (from [UniProt entry Q12809](https://www.uniprot.org/uniprotkb/Q12809))
- [KMT2A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=KMT2A) - Histone-lysine N-methyltransferase 2A (from [UniProt entry Q03164](https://www.uniprot.org/uniprotkb/Q03164))
- [LMNA](https://www.genecards.org/cgi-bin/carddisp.pl?gene=LMNA) - Prelamin-A/C (from [UniProt entry P02545](https://www.uniprot.org/uniprotkb/P02545))
- [LTA4H](https://www.genecards.org/cgi-bin/carddisp.pl?gene=LTA4H) - Leukotriene A-4 hydrolase (from [UniProt entry P09960](https://www.uniprot.org/uniprotkb/P09960))
- [MAOB](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MAOB) - Amine oxidase [flavin-containing] B (from [UniProt entry P27338](https://www.uniprot.org/uniprotkb/P27338))
- [MAPK1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MAPK1) - Mitogen-activated protein kinase 1 (from [UniProt entry P28482](https://www.uniprot.org/uniprotkb/P28482))
- [MDH1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MDH1) - Malate dehydrogenase, cytoplasmic (from [UniProt entry P40925](https://www.uniprot.org/uniprotkb/P40925))
- [MMP2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MMP2) - 72 kDa type IV collagenase (from [UniProt entry P08253](https://www.uniprot.org/uniprotkb/P08253))
- [MMP3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MMP3) - Stromelysin-1 (from [UniProt entry P08254](https://www.uniprot.org/uniprotkb/P08254))
- [MMP9](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MMP9) - Matrix metalloproteinase-9 (from [UniProt entry P14780](https://www.uniprot.org/uniprotkb/P14780))
- [MTCA1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MTCA1) - Beta-carbonic anhydrase 1 (from [UniProt entry P9WPJ7](https://www.uniprot.org/uniprotkb/P9WPJ7))
- [MTCA2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=MTCA2) - Carbonic anhydrase 2 (from [UniProt entry P9WPJ9](https://www.uniprot.org/uniprotkb/P9WPJ9))
- [NCE103](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NCE103) - Carbonic anhydrase (from [UniProt entry Q5AJ71](https://www.uniprot.org/uniprotkb/Q5AJ71))
- [NFE2L2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NFE2L2) - Nuclear factor erythroid 2-related factor 2 (from [UniProt entry Q16236](https://www.uniprot.org/uniprotkb/Q16236))
- [NFKB1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NFKB1) - Nuclear factor NF-kappa-B p105 subunit (from [UniProt entry P19838](https://www.uniprot.org/uniprotkb/P19838))
- [NFO](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NFO) - Probable endonuclease 4 (from [UniProt entry A0AIQ1](https://www.uniprot.org/uniprotkb/A0AIQ1))
- [NPSR1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NPSR1) - Neuropeptide S receptor (from [UniProt entry Q6W5P4](https://www.uniprot.org/uniprotkb/Q6W5P4))
- [NS1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NS1) - Influenza virus NS1A-binding protein (from [UniProt entry Q9Y6Y0](https://www.uniprot.org/uniprotkb/Q9Y6Y0))
- [P2RX4](https://www.genecards.org/cgi-bin/carddisp.pl?gene=P2RX4) - P2X purinoceptor 4 (from [UniProt entry Q99571](https://www.uniprot.org/uniprotkb/Q99571))
- [P2RX7](https://www.genecards.org/cgi-bin/carddisp.pl?gene=P2RX7) - P2X purinoceptor 7 (from [UniProt entry Q99572](https://www.uniprot.org/uniprotkb/Q99572))
- [PDE3A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=PDE3A) - cGMP-inhibited 3',5'-cyclic phosphodiesterase 3A (from [UniProt entry Q14432](https://www.uniprot.org/uniprotkb/Q14432))
- [PFK](https://www.genecards.org/cgi-bin/carddisp.pl?gene=PFK) - ATP-dependent 6-phosphofructokinase (from [UniProt entry P00512](https://www.uniprot.org/uniprotkb/P00512))
- [PMP22](https://www.genecards.org/cgi-bin/carddisp.pl?gene=PMP22) - Peroxisomal membrane protein 2 (from [UniProt entry Q9NR77](https://www.uniprot.org/uniprotkb/Q9NR77))
- [PTH1R](https://www.genecards.org/cgi-bin/carddisp.pl?gene=PTH1R) - Parathyroid hormone/parathyroid hormone-related peptide receptor (from [UniProt entry Q03431](https://www.uniprot.org/uniprotkb/Q03431))
- [REP](https://www.genecards.org/cgi-bin/carddisp.pl?gene=REP) - Replicase polyprotein 1ab (from [UniProt entry P0C6W7](https://www.uniprot.org/uniprotkb/P0C6W7))
- [RGS4](https://www.genecards.org/cgi-bin/carddisp.pl?gene=RGS4) - Regulator of G-protein signaling 4 (from [UniProt entry P49798](https://www.uniprot.org/uniprotkb/P49798))
- [RORC](https://www.genecards.org/cgi-bin/carddisp.pl?gene=RORC) - Nuclear receptor ROR-gamma (from [UniProt entry P51449](https://www.uniprot.org/uniprotkb/P51449))
- [SAMHD1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SAMHD1) - Deoxynucleoside triphosphate triphosphohydrolase SAMHD1 (from [UniProt entry Q9Y3Z3](https://www.uniprot.org/uniprotkb/Q9Y3Z3))
- [SCN10A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SCN10A) - Sodium channel protein type 10 subunit alpha (from [UniProt entry Q9Y5Y9](https://www.uniprot.org/uniprotkb/Q9Y5Y9))
- [SCN1A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SCN1A) - Sodium channel protein type 1 subunit alpha (from [UniProt entry P35498](https://www.uniprot.org/uniprotkb/P35498))
- [SCN2A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SCN2A) - Sodium channel protein type 2 subunit alpha (from [UniProt entry Q99250](https://www.uniprot.org/uniprotkb/Q99250))
- [SCN3A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SCN3A) - Sodium channel protein type 3 subunit alpha (from [UniProt entry Q9NY46](https://www.uniprot.org/uniprotkb/Q9NY46))
- [SCN4A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SCN4A) - Sodium channel protein type 4 subunit alpha (from [UniProt entry P35499](https://www.uniprot.org/uniprotkb/P35499))
- [SCN5A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SCN5A) - Sodium channel protein type 5 subunit alpha (from [UniProt entry Q14524](https://www.uniprot.org/uniprotkb/Q14524))
- [SCN8A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SCN8A) - Sodium channel protein type 8 subunit alpha (from [UniProt entry Q9UQD0](https://www.uniprot.org/uniprotkb/Q9UQD0))
- [SCN9A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SCN9A) - Sodium channel protein type 9 subunit alpha (from [UniProt entry Q15858](https://www.uniprot.org/uniprotkb/Q15858))
- [SLC22A1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC22A1) - Solute carrier family 22 member 1 (from [UniProt entry O15245](https://www.uniprot.org/uniprotkb/O15245))
- [SLC22A6](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC22A6) - Solute carrier family 22 member 6 (from [UniProt entry Q4U2R8](https://www.uniprot.org/uniprotkb/Q4U2R8))
- [SLC6A3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC6A3) - Sodium-dependent dopamine transporter (from [UniProt entry Q01959](https://www.uniprot.org/uniprotkb/Q01959))
- [SLCO1B1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLCO1B1) - Solute carrier organic anion transporter family member 1B1 (from [UniProt entry Q9Y6L6](https://www.uniprot.org/uniprotkb/Q9Y6L6))
- [SLCO1B3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLCO1B3) - Solute carrier organic anion transporter family member 1B3 (from [UniProt entry Q9NPD5](https://www.uniprot.org/uniprotkb/Q9NPD5))
- [SLCO2B1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLCO2B1) - Solute carrier organic anion transporter family member 2B1 (from [UniProt entry O94956](https://www.uniprot.org/uniprotkb/O94956))
- [SMN2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SMN2) - Survival motor neuron protein (from [UniProt entry Q16637](https://www.uniprot.org/uniprotkb/Q16637))
- [TDP1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TDP1) - Tyrosyl-DNA phosphodiesterase 1 (from [UniProt entry Q9NUW8](https://www.uniprot.org/uniprotkb/Q9NUW8))
- [THRB](https://www.genecards.org/cgi-bin/carddisp.pl?gene=THRB) - Thyroid hormone receptor beta (from [UniProt entry P10828](https://www.uniprot.org/uniprotkb/P10828))
- [TP53](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TP53) - Cellular tumor antigen p53 (from [UniProt entry P04637](https://www.uniprot.org/uniprotkb/P04637))
- [TSHR](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TSHR) - Thyrotropin receptor (from [UniProt entry P16473](https://www.uniprot.org/uniprotkb/P16473))
- [TXNRD1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TXNRD1) - Thioredoxin reductase 1, cytoplasmic (from [UniProt entry Q16881](https://www.uniprot.org/uniprotkb/Q16881))
- [UGT1A9](https://www.genecards.org/cgi-bin/carddisp.pl?gene=UGT1A9) - UDP-glucuronosyltransferase 1A9 (from [UniProt entry O60656](https://www.uniprot.org/uniprotkb/O60656))
- [UGT2B15](https://www.genecards.org/cgi-bin/carddisp.pl?gene=UGT2B15) - UDP-glucuronosyltransferase 2B15 (from [UniProt entry P54855](https://www.uniprot.org/uniprotkb/P54855))
- [VPR](https://www.genecards.org/cgi-bin/carddisp.pl?gene=VPR) - Minor extracellular protease Vpr (from [UniProt entry P29141](https://www.uniprot.org/uniprotkb/P29141))
