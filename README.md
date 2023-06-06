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
$ conda create -c rdkit -n rdkit-env rdkit
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

This should be a trivial process. Just copy the `csanno.py` files to a new folder in your woring directory: 

```
$ mkdir qsartools
$ cp [source code path]/*.py qsartools/.
$ cd qsartools
```

### SAR files

SAR files are the required format for most of the operations. SAR files are text files that should store activity data for a given target for a set of molecules and each line represents one molecule. Columns are separated by `<tab>`. The basic structur of a SAR file is

* The first column is an alphanumeric identifier for the molecule (should be unique);
* The second column contains the activity registered for that molecule (may be binary qualitative)
* The third column contains the molecule SMILES

Here is a sample SAR file, where the identifiers are the ChEMBL IDs

```
CHEMBL158973	0.5000	CN(C)CCSC(C)(C)C
CHEMBL476516	0.1880	NCCc1ccc(Br)cc1
CHEMBL309689	0.0000	Oc1noc2c1CCNCC2
...
CHEMBL2331792	0.0000	COc1noc2c1CCNCC2
CHEMBL2331804	0.0000	Cn1oc2c(c1=O)CCNCC2
CHEMBL2436555	0.2200	NC1=Nc2ccc(Cl)cc2CN1
```

## Asking for Help

`csanno` has a built in helpthat clarifies all its options

```sh
$ python csanno.py -help
```
producing
```
CSANNO - (C) 2019/2023 - Andre O. Falcao DI/FCUL version 0.3.20230605
Usage: This is a Python3 tool and requires an enviroment where RDkit and requests are installed.
       To run type python csanno.py  -in [input .sar file] [options]

Control parameters:
    -in file_name - the data set to annotate (required) (.sar format)
    -sim thr - similarity threshold [if absent (default) no similarity search will be performed]
    -mol [molecule in SMILES format] - molecule for which analysis is going to be performed.
         Incompatible with -in option
    -report [AD]: A - aggregate report;
                  D - Detailed report (one line for each molecule)
                  Both options simultaneusly are permitted and both reports are produced
    -search [ANU] A - searches for actives in the database (default);
                  N - searches non-actives;
                  U - searches unknowns (compounds for which activity is not determined)
                  All options can be included simultaneously
    -ofield [GUO] G - outputs the unique gene (default);
                  U - outputs the uniprot id;
                  O - outputs the organism
                  All options can be included simultaneously

Output Control Options:
    -joint_aggs - produces a report with conjoined similarities and exact matches
    -silent - no intermediate output at all
    -help - this screen
    -nofiles - no files are created
    -toscreen - writes output to screen
    -pickle - writes a Python pickle for easy posterior analysis (see internal documentation)
```

Essentially it allows to search for ChEMBL structural similars and ellaborating a report out of it


## Molecular activity annotation

`csanno` is a molecular annotator that uses ChEMBL to search to identify plausible targets. The question we are trying to solve in this program could be better defined with a simple example. Imagine we have one molecule and we want to identify probable targets in which it may be active. There are two possible situations: 

1. The molecule is present in ChEMBL and therefore all the targets for which it was found active are retrieved
2. The molecule is not found in ChEMBL. And in that case what it will do is that it will look for similar molecules (using ChEMBL's similarity search) and will output the found targets and the number of molecule similars that were actually found active there

The tool also has a input file, option where for instance we might group a set of molecules with similar physiological activity, to try to identify common or most likely targets. All of this would be actually possible to do manually using ChEMBL online, but it would be tedious and extremely prone to errors. 


### Simple single molecule annotation

The simplest way of running the csanno is by testing only one molecule that can function as an input

Using the `-mol` option writes the result directly to the screen. Therefore running the app with zolpidem (SMILES:CN(C)C(=O)Cc1c(nc2ccc(C)cn12)c3ccc(C)cc3)

```sh
$ python csanno.py -mol CN(C)C(=O)Cc1c(nc2ccc(C)cn12)c3ccc(C)cc3
```

will produce the following output:

```
TSPO                 1  1.0000
SLCO1B3              1  1.0000
SLCO1B1              1  1.0000
LMNA                 1  1.0000
GABRA3               1  1.0000
GABRA2               1  1.0000
GABRA1               1  1.0000
```

This result is a table which must be read this way. In the first column it is the gene symbol associated with a given assay for which the molecule has been proven active. It is important to notice that **activity** here means that a positive response has been found for it in at least one assay, whether it is Activation, Ki, IC50 or simple inhibition. If a positive result has been found, then it is reported. The next column refers to how many molecules were actually measured as active in that target. As we have only used one molecule, the only possible result is 1. Finally, the last column with values of 1.0000, refers to the ratio of the input molecules for which there were found actives. Again as we have tested only one molecule, this is the only possible output

### Simple multiple molecule annotation

To test for the multiple targets we have created a simple file with 5 molecules, all known anti depressants to show the common binding profiles 

```
Sertraline	NA	ClC1=CC=C([C@H]2C3=C([C@H](CC2)NC)C=CC=C3)C=C1Cl
Citalopran	NA	Fc1ccc(cc1)C3(OCc2cc(C#N)ccc23)CCCN(C)C
Fluoxetine	NA	CNCCC(c1ccccc1)Oc2ccc(cc2)C(F)(F)F
Tradozone	NA	Clc4cccc(N3CCN(CCCN1/N=C2/C=C\C=C/N2C1=O)CC3)c4
Venlafaxine	NA	OC2(C(c1ccc(OC)cc1)CN(C)C)CCCCC2
```

Running the following command

```sh
$ python csanno.py -in anti-dep.sar -sim 0.7 -report AD 
```

will produce 2 files. The name of those files is created based on the .sar file. So using `anti-dep.sar` will produce:

* `anti-dep_anno_count_T0.txt` - This is a Tier 0 file and will contain the absolute and relative frequencies of the molecules in the input file that appear to be active in specific targets
* `anti-dep_anno_count_T1.txt` - This is a Tier 1 file and will contain the the absolute and relative frequencies of the molecules **that are similar** to the molecules in the input file to each specific target in which they were found active

The tier 0 file thus refers to the known targets in single-gene essays that are know to bind to these proteins. `anti-dep_anno_count_T0.txt` will look something like this:


```
SLC6A4	    5	1.0000
SLC6A3	    4	0.8000
SLC6A2	    4	0.8000
SIGMAR1    4	0.8000
HTR2C	    4	0.8000
REP	    3	0.6000
NET	    3	0.6000
KCNH2	    3	0.6000
HTR2B	    3	0.6000
```

This file is actually 52 lines long but only the top 8 are shown. It shows which are the individual targets common in all of them. As we can see, [SLC6A4](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC6A4) (the serotonin transporter) is common to all 5 molecules, followed by [SLC6A3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC6A3) and [SLC6A2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC6A2) the dopamine and noradrenaline transporter, respectively. Also the the [Sigma1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SIGMA1) Receptor and the [Serotonin 2C receptor](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR2C) are known to be affected by 4 out of 5 molecules

On the other hand we can look at `anti-dep_anno_count_T1.txt` for potential targets, as these refer to the activity profiles of similar molecules (in this case 70% similar) 

```
SLC6A4	   64	0.6737
SLC6A2	   29	0.3053
SLC6A3	   24	0.2526
REP	    8	0.0842
ADRA1A	    8	0.0842
NET	    7	0.0737
HTR2A	    6	0.0632
SIGMAR1    4	0.0421
```
of the 183 molecules found, the top spots were essentially the same, but the [ADRA1A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ADRA1A) Adrenergic A1 receptor and the [NET](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NET), the norepinephrine transporter, appear with a high representativity. Troubling might be the presence of [HTR2A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR2A), in 6% of all matches. 

To have a better look at which of the base molecules are associated, we may run the above command again with the following option `-report AD` This will create two more files that stand for the individual targets each molecule was active on. The format is a standard 'transaction type' file with one line for each molecule, followed by the found targets, named similarly to the other two files created:

* `anti-dep_annotations_T0.txt` - binding targets to each molecule in the input file
* `anti-dep_annotations_T1.txt` - binding targets to each molecule similar to the molecules in the input file.

Therefore `anti-dep_annotations_T0.txt` with the above options will look like this (end of lines arranged for this report, in the original, only one line per molecule is produced):

```
Sertraline:	ADRA2B	CACNA1C	KMT2A	SIGMAR1	CYP3A4	ADRA2A	HTR2B	REP	CYP2C19	
              CHRM1	KCNH2	CHRM4	ADRA2C	MC5R	CHRM5	HTR2A	SLC6A3	HTR2C	CYP2D6	ADRA1B	
              TMEM97	SLC6A2	CHRM2	SLC6A4
Citalopran:	TMEM97	ADRA1B	SLC6A2	HRH1	KCNH2	ADRA1A	NET	SIGMAR1	SMN2	ADRA1D	
              SLC6A3	HTR2C	REP	HTR2B	SLC6A4	SLC22A1
Fluoxetine:	ABCB1	ADRA2B	CACNA1C	KCNK9	SIGMAR1	NET	HRH3	ACHE	CYP3A4	
              ADRA2A	SLCO1B1	HTR3A	CYP2C19	CHRM1	CYP2C9	HRH1	KCNH2	LMNA	
              KCNK2	CHRM5	CYP2B6	HTR6	CHRM3	HTR2A	DRD2	SLC6A3	SLCO2B1	HTR2C	
              CYP1A2	CYP2D6	CYP2C8	SLC6A2	SLCO1B3	SLC6A4
Tradozone:	ADRA2C	HTR1B	ADRA2B	SLC6A4	HRH1	ADRA1A	SIGMAR1	HTR2A	DRD3	ALB	
              ADRA1D	DRD2	FAAH	HTR2C	ADRA2A	HTR2B	HTR1A	ADRA1B
Venlafaxine:	SLC6A2	ADRA1A	NET	LMNA	SLC6A3	REP	SLC6A4
```

Similarly `anti-dep_annotations_T1.txt` with the above options will look like this:


```
Sertraline:	NFKB1	TMEM97	SLCO1B3	SLCO1B1	REP	HTR2A	KCNH2	CYP3A4	ADRA1B	
              CHRM2	SLC6A2	NET	NS1	MAPK1	SLC6A3	ADRA2A	SIGMAR1	CYP1A2	FFP	
              ADRA2C	KMT2A	TP53	MC5R	CACNA1C	HTR2B	CYP2D6	CHRM4	SLC6A4	ADRA2B	
              CHRM1	CHRM5	HTR2C	CYP2C19
Citalopran:	NFKB1	TMEM97	SLCO1B3	SLCO1B1	REP	KCNH2	ADRA1B	CYP3A4	SLC6A2	
              SLC22A1	ADRA1D	ADRA1A	NET	MAPK1	SLC6A3	SIGMAR1	ALOX15	SMN2	
              CYP2C9	HTR2B	SLC6A4	CHRM1	HRH1	LMNA	HTR2C	CYP2C19
Fluoxetine:	CYP2B6	NFKB1	HTR3A	CYP2C8	SLCO1B3	SLCO1B1	HTR2A	KCNK9	KCNH2	
              DRD2	CYP3A4	KCNK2	SLC6A2	TSHR	NET	ADRA2A	SLC6A3	SIGMAR1	SLCO2B1
              CYP1A2	HTR6	HRH3	SMN2	CYP2C9	CACNA1C	BLM	CYP2D6	ACHE	ABCB1	
              RORC	SLC6A4	ADRA2B	CHRM1	HRH1	LMNA	CHRM3	CHRM5	HTR2C	AMPC	HPGD	
              CYP2C19
Tradozone:	FAAH	SLCO1B3	SLCO1B1	REP	HTR2A	DRD2	ADRA1B	SLC6A2	ADRA1D	
              ADRA1A	HTR1A	ADRA2A	SIGMAR1	ADRA2C	HTR2B	ADRA2B	SLC6A4	HRH1	HTR2C	HTR1B	ALB	DRD3	HLA-A
Venlafaxine:	SLC6A4	SLC6A2	LMNA	ADRA1A	NET	SLC6A3	REP
```

Please note that these results are the targets **that do not apear on Tier 0**, meaning possible targets never before tested

As an example of a possible conclusion is that in some of the 70% similar molecules to Sertraline some of them were verified active on NET and CYP2C9, strongly suggesting that Sertraline itself, although not registered as active in these two targets, might actually be so. Other types of observations might be held for other molecules.
