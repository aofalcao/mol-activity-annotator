# CSANNO: a molecular activity annotator

`csanno` is a molecular activity-evidence annotator that uses ChEMBL to retrieve experimental records for a query molecule and, optionally, for structurally similar molecules. Two separate evidence tiers are available:

1. **Tier 0 (direct evidence):** activity records retrieved for the input molecule through exact-structure identification.
2. **Tier 1 (analogue-derived evidence):** activity records retrieved for structurally similar ChEMBL molecules when similarity searching is requested with `-sim`. This can be performed whether or not an exact match was found in Tier 0.

CSANNO distinguishes active, non-active, and uncertain evidence according to configurable rules. Tier 1 evidence describes measured activity of neighbouring molecules, not direct measurements for the query molecule. The tool also accepts input files containing sets of molecules, allowing common evidence patterns to be examined across a pharmacological or chemical class. These operations could be performed manually using ChEMBL, but would be tedious and error-prone.

`csanno` currently relies on [ChEMBL web services](https://pubmed.ncbi.nlm.nih.gov/28602100/), so a working internet connection is needed for remote queries. A ChEMBL-compatible server may also be specified with `-chembl\_server`.

The annotation profile of a molecule is a report of the experimental evidence retrieved for that structure (Tier 0) and, if requested, for structurally similar molecules (Tier 1), retaining the distinction between the two sources of evidence.

## CSANNO requirements

`CSANNO` is Python software that relies on the following external libraries (also listed in `requirements.txt`):

* requests
* rdkit
* PyYAML

Of these, RDKit is often the trickiest to install. A conda environment is recommended, following the [official RDKit installation instructions](https://www.rdkit.org/docs/Install.html).

First create an environment `rdkit-env`:

```sh
$ conda create -c conda-forge -n rdkit-env rdkit
```

Then activate it:

```sh
$ conda activate rdkit-env
```

### Using `requirements.txt`

The file `requirements.txt` lists the required dependencies. Install them in a suitable Python environment with:

```sh
$ python3 -m pip install -r requirements.txt
```

It is recommended that RDKit be installed in a dedicated Python environment.

### Installing CSANNO

Clone the repository or copy the `.py` and `.yaml` files to a working directory. The default YAML activity-rule profile must be available to the program.

```sh
$ git clone https://github.com/aofalcao/mol-activity-annotator.git
$ cd mol-activity-annotator
```

### Molecular standardization (`-standardize`)

Molecular standardization is **optional and disabled by default** to preserve the original CSANNO behaviour. When `-standardize` is specified, each input SMILES is processed with RDKit's `MolStandardize` operations in the following order: `Cleanup`, `FragmentParent`, `Uncharger`, and `IsotopeParent`. A canonical isomeric SMILES is generated from the processed structure.

This same processed structure is used to generate the InChIKey for **Tier 0** exact lookup and as the SMILES submitted for **Tier 1** similarity retrieval. Without `-standardize`, the input structure is used without these additional standardization operations (although RDKit still parses it and generates its InChI/InChIKey for exact lookup). Standardization can remove counterions, select a parent fragment, neutralize charges, or remove isotope labels, and may therefore change retrieved matches. It does not perform a comprehensive tautomer normalization or guarantee identical processing to ChEMBL's internal standardization.

For example, append `-standardize` to a normal command when this processing is desired:

```sh
$ python csanno.py -mol 'CC(=O)\[O-].\[Na+]' -sim 0.7 -standardize
```


### SAR files

SAR files are the required format for file-based operations. They are tab-delimited text files with one molecule per line:

* The first column is a unique alphanumeric identifier for the molecule.
* The second column contains a user-supplied activity or class value (which is retained for traceability but does not determine the database-derived evidence classification).
* The third column contains the molecule SMILES.

Here is a sample SAR file containing antidepressants. The middle column is `NA`, meaning that the user-supplied activity is unknown or not relevant:

```
Sertraline	NA	ClC1=CC=C(\[C@H]2C3=C(\[C@H](CC2)NC)C=CC=C3)C=C1Cl
Citalopran	NA	Fc1ccc(cc1)C3(OCc2cc(C#N)ccc23)CCCN(C)C
Fluoxetine	NA	CNCCC(c1ccccc1)Oc2ccc(cc2)C(F)(F)F
Tradozone	NA	Clc4cccc(N3CCN(CCCN1/N=C2/C=C\\C=C/N2C1=O)CC3)c4
Venlafaxine	NA	OC2(C(c1ccc(OC)cc1)CN(C)C)CCCCC2
```

## Asking for Help

`csanno` has built-in help describing its options:

```sh
$ python csanno.py -help
```

producing:

```
usage: csanno.py \[-in FILE | -mol SMILES] \[-sim THRESHOLD] \[-search ANU] \[-report {A,D,AD}] \[-ofield GUO]
                 \[-chembl\_server URL] \[-rules YAML\_FILE] \[-standardize] \[-out NAME] \[-outfmt {markdown,md,text,txt}]
                 \[-legacy\_reports] \[-joint\_aggs] \[-glossary] \[-digest] \[-silent] \[-nofiles] \[-toscreen] \[-pickle]
                 \[-help] \[--version]
Annotate molecules with direct ChEMBL evidence (Tier 0) and, optionally, evidence from structurally similar ChEMBL
molecules (Tier 1).
options:
  -in FILE              input data set in three-column tab-separated .sar format
  -mol SMILES           analyse one molecule supplied as a SMILES string
  -sim THRESHOLD        enable Tier 1 similarity searching at this threshold
  -search ANU           evidence classes to include: A, N and/or U (default: A)
  -report {A,D,AD}      report sections: A, D or AD (default: AD)
  -ofield GUO           target fields to output: G, U and/or O (default: G)
  -chembl\_server URL    ChEMBL-compatible data server
  -rules YAML\_FILE      YAML file defining activity classification rules
  -standardize          standardize input molecules before Tier 0 and Tier 1 searches
  -out NAME             output file root
  -outfmt {markdown,md,text,txt}
                        readable report format (default: markdown)
  -legacy\_reports, -legacy
                        also write the historical split report files
  -joint\_aggs           with -sim and legacy reports, write joined T0+T1 reports
  -glossary             include a gene glossary
  -digest               include the report digest
  -silent               suppress intermediate output
  -nofiles              do not write output files
  -toscreen             write the readable report to the screen
  -pickle               write a Python pickle for later analysis
  -help, --help         show this help message and exit
  --version             show program's version number and exit
Output Control Options:
    By default, CSANNO writes one readable Markdown report (<root>\_report.md) and one
    JSON report (<root>\_results.json) containing metadata and raw data.
    -out file\_name - output file name
    -outfmt \[markdown|text] - writes the readable report as .md (default) or .txt
    -legacy\_reports - also writes the old split T0/T1 aggregate/detail text files
    -joint\_aggs - with -sim and -legacy\_reports, also produces joined T0+T1 legacy files
    -glossary - includes a glossary of the targets found with links
    -digest - includes an abridged analysis of Tier 0 and Tier 1 reports
    -silent - no intermediate output at all
    -help - this screen
    -nofiles - no files are created
    -toscreen - writes the readable report to screen
    -pickle - writes a Python pickle for later analysis (see internal documentation)
```

## Molecular activity annotation

CSANNO retrieves ChEMBL activity records associated with exact query structures (Tier 0), and optionally with structurally similar compounds (Tier 1). These reports can be generated for individual molecules or for sets of molecules. The second tier is enabled explicitly with `-sim`, irrespective of whether an exact ChEMBL match was found. Results are organized by target label and evidence category; an annotation is evidence retrieved under selected rules, not a new prediction of target activity.

### Simple single molecule annotation

The simplest way of running CSANNO is to supply one molecule. For example, with zolpidem (SMILES: `CN(C)C(=O)Cc1c(nc2ccc(C)cn12)c3ccc(C)cc3`):

```sh
$ python csanno.py -mol 'CN(C)C(=O)Cc1c(nc2ccc(C)cn12)c3ccc(C)cc3' -sim 0.7 -out data/zolpidem-A
```

This produces a report like the [zolpidem example](https://github.com/aofalcao/mol-activity-annotator/blob/main/Reports/zolpidem-A_report.md).

### Simple multiple molecule annotation

To analyse several molecules, use the sample five-antidepressant SAR file above. For example:

```sh
$ python csanno.py -in data/anti-dep.sar -sim 0.7 -report AD -digest -glossary
```

This command:

* Analyses the five query molecules using exact-match (Tier 0) evidence.
* Retrieves similar ChEMBL molecules at the selected threshold of 70% (`-sim 0.7`) and reports their evidence separately (Tier 1).
* Creates aggregate and detailed reports (`-report AD`), including counts of input molecules or neighbours supporting each target label.
* Includes a concise report digest (`-digest`).
* Includes a glossary of target names and links (`-glossary`).

To apply the optional standardization procedure to this example, add `-standardize` to the command. The [example report](https://github.com/aofalcao/mol-activity-annotator/blob/main/Reports/anti-dep_report.md) illustrates the reporting format; its precise outputs may depend on the ChEMBL data available at query time.

