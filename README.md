# mol-activity-annotator

`csanno` is a molecular annotator that uses ChEMBL to search to identify plausible targets. The question we are trying to solve in this program could be better defined with a simple example. Imagine we have one molecule and we want to identify probable targets in which it may be active. There are two possible situations: 

1. The molecule is present in ChEMBL and therefore all the targets for which it was found active are retrieved
2. The molecule is not found in ChEMBL. And in that case what it will do is that it will look for similar molecules (using ChEMBL's similarity search) and will output the found targets and the number of molecule similars that were actually found active there

The tool also has a input file, option where for instance we might group a set of molecules with similar physiological activity, to try to identify common or most likely targets. All of this would be actually possible to do manually using ChEMBL online, but it would be tedious and extremely prone to errors. 


`csanno` is heavily reliant on [ChEMBL web services](https://pubmed.ncbi.nlm.nih.gov/28602100/), therefore it is fundamental to have a good internet connection. Alternatively if a local ChEMBL installation is available it can be easily used with increased performance benefits.

The annotation profile of a molecule is essentially a report of the evidence that ChEMBL has for that specific structure or, if required, of the molecules that are very similar to it.


## CSANNO requirements

`CSANNO` is essentially Python software that relies on a number of external libraries. Namely:

* numpy
* requests
* rdkit

Of those, rdkit is the trickiest to install. Anaconda makes the rdkit installation a breeze and this is the process recommended, although it might be possible to use a [docker container](https://hub.docker.com/search?q=rdkit&type=image) with rdkit already installed. To install rdkit with anaconda we can follow the [procedure detailed in the official rdkit distribution](https://www.rdkit.org/docs/Install.html):

First we create an environment `rdkit-env` like this:

```sh
$ conda create -c conda-forge -n rdkit-env rdkit
```

and then to activate we just activate the environment to start working

```sh
$ conda activate rdkit-env
```

Windows users can omit the `conda` and just use `activate rdkit-env` at the windows command line prompt


With Anaconda, numpy is installed by default and installing the others can be performed simply with: 

```sh
$ conda install requests
```

Just make sure these are installed within the `rdkit` environment. 

### Installing CSANNO

This should be a trivial process. Either use git or just copy the `.py` files to a new folder in your woring directory: 

```
$ mkdir csanno
$ cp [source code path]/*.py csanno/.
$ cd csanno
```

### SAR files

SAR files are the required format for file based operations. SAR files are text files that should store activity data for a given target for a set of molecules and each line represents one molecule. Columns are separated by `<tab>`. The basic structur of a SAR file is

* The first column is an alphanumeric identifier for the molecule (should be unique);
* The second column contains the activity registered for that molecule (may be binary qualitative)
* The third column contains the molecule SMILES

Here is a sample SAR file, where the identifiers are the Antidepressants. The middle column is NA, meaning that the activity is unknown or not relevant

```
Sertraline	NA	ClC1=CC=C([C@H]2C3=C([C@H](CC2)NC)C=CC=C3)C=C1Cl
Citalopran	NA	Fc1ccc(cc1)C3(OCc2cc(C#N)ccc23)CCCN(C)C
Fluoxetine	NA	CNCCC(c1ccccc1)Oc2ccc(cc2)C(F)(F)F
Tradozone	NA	Clc4cccc(N3CCN(CCCN1/N=C2/C=C\C=C/N2C1=O)CC3)c4
Venlafaxine	NA	OC2(C(c1ccc(OC)cc1)CN(C)C)CCCCC2
```

## Asking for Help

`csanno` has a built in helpthat clarifies all its options

```sh
$ python csanno.py -help
```
producing
```
(rdkit-env) C:\Users\aofal\Desktop\Projectos\CSANNO>python csanno.py -help
CSANNO - (C) 2019/2026 - Andre O. Falcao BioISI - DI/FCUL version 0.5.20260515
Usage: This is a Python tool and requires an enviroment where RDkit and requests are installed.
       To run type: python csanno.py  -in [input .sar file] [options] OR
                    python csanno.py  -mol [molecule SMILES] [options]
Control parameters:
    -in file_name - the data set to annotate (required) (.sar format)
    -sim thr - similarity threshold [if absent (default) no similarity search will be performed]
    -mol [molecule in SMILES format] - molecule for which analysis is going to be performed.
         Incompatible with -in option
    -report [AD]: A - aggregate report;
                  D - Detailed report (one line for each molecule)
                  Both options simultaneously are permitted and both reports are produced
    -search [ANU] A - searches for actives in the database (default);
                  N - searches non-actives in exact ChEMBL matches;
                  U - searches unknowns in exact ChEMBL matches;
                  Similarity reports use the same search option
                  All options can be included simultaneously
    -ofield [GUO] G - outputs the unique gene (default);
                  U - outputs the uniprot id;
                  O - outputs the organism;
                  All options can be included simultaneously
    -chembl_server URL - ChEMBL-compatible data server [default: https://www.ebi.ac.uk]
    -rules YAML_file - YAML file defining active/non-active/uncertain classification rules

Output Control Options:
    By default, CSANNO writes one readable Markdown report (<root>_report.md) and one
    JSON report (<root>_results.json) containing metadata and raw data.
    -out file_name - output file name
    -outfmt [markdown|text] - writes the readable report as .md (default) or .txt
    -legacy_reports - also writes the old split T0/T1 aggregate/detail text files
    -joint_aggs - with -sim and -legacy_reports, also produces joined T0+T1 legacy files
    -glossary - includes a glossary of the targets found with links
    -digest - includes an abbridged analysis of Tier 0 and Tier 1 reports
    -silent - no intermediate output at all
    -help - this screen
    -nofiles - no files are created
    -toscreen - writes the readable report to screen
    -pickle - writes a Python pickle for easy posterior analysis (see internal documentation)
```

Essentially it allows to search for structural similars (Currently using ChEMBLE REST APIs only) and ellaborating a report out of it


## Molecular activity annotation

`csanno` is a molecular annotator that uses A molecular  activity database such as ChEMBL to search to identify plausible targets. The question we are trying to solve in this program could be better defined with a simple example. Imagine we have one molecule and we want to identify probable targets in which it may be active. There are two possible situations: 

1. The molecule is present in the database and therefore all the targets for which it was found active are retrieved
2. The molecule is not found in the database. And in that case what it will do is that it will look for similar molecules (using ChEMBL's similarity search) and will output the found targets and the number of molecule similars that were actually found active there

The tool also has a input file, option where for instance we might group a set of molecules with similar physiological activity, to try to identify common or most likely targets. All of this would be actually possible to do manually browsing and querying the database, but it would be tedious and extremely prone to errors. 


### Simple single molecule annotation

The simplest way of running the csanno is by testing only one molecule that can function as an input

Therefore running the app with zolpidem (SMILES:CN(C)C(=O)Cc1c(nc2ccc(C)cn12)c3ccc(C)cc3)

```sh
$ python csanno.py -mol CN(C)C(=O)Cc1c(nc2ccc(C)cn12)c3ccc(C)cc3 -sim 0.7 -out data\zolpidem-A
```

will produce the [following output](https://github.com/aofalcao/mol-activity-annotator/blob/main/Reports/zolpidem-A_report.md)


### Simple multiple molecule annotation

To test for the multiple targets we will use the above mentioned file 5 molecules, all known anti depressants to show the common binding profiles 

```
Sertraline	NA	ClC1=CC=C([C@H]2C3=C([C@H](CC2)NC)C=CC=C3)C=C1Cl
Citalopran	NA	Fc1ccc(cc1)C3(OCc2cc(C#N)ccc23)CCCN(C)C
Fluoxetine	NA	CNCCC(c1ccccc1)Oc2ccc(cc2)C(F)(F)F
Tradozone	NA	Clc4cccc(N3CCN(CCCN1/N=C2/C=C\C=C/N2C1=O)CC3)c4
Venlafaxine	NA	OC2(C(c1ccc(OC)cc1)CN(C)C)CCCCC2
```

```sh
$ python csanno.py -in data\anti-dep.sar -sim 0.7 -report AD -digest -glossary
```
The above basically says:
* Analyze the file `data\anti-dep.sar` by the direct evidence there is for the whole molecules
* Check all molecules in the database that are currently similar to it at 70% (`-sim 0.7`)
* Create an aggregate and detailed report (`-report AD`) getting, for each molecule and each tier (direct vs analogue evidence), the actual targets found, relative to the targets, how many analogues were identified there
* create a digest - interpretable analysis of the results (`-digest`)
* create a glossary of the targets found with links to Chembl and Uniprot (`-glossary`)
  

Results can be [found here](https://github.com/aofalcao/mol-activity-annotator/blob/main/Reports/anti-dep_report.md)



