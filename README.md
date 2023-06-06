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
$ python csanno.py -in anti-dep.sar -sim 0.8 
```

will produce 2 files. The name of those files is created based on the .sar file. So using `anti-dep.sar` will produce:

* `anti-dep_anno_count_T0.txt` - This is a Tier 0 file and will contain the absolute and relative frequencies of the molecules in the input file that appear to be active in specific targets
* `anti-dep_anno_count_T1.txt` - This is a Tier 1 file and will contain the the absolute and relative frequencies of the molecules **that are similar** to the molecules in the input file to each specific target in which they wer found active

The tier 0 file thus refers to the known targets in single-gene essays that are know to bind to these proteins. `anti-dep_anno_count_T0.txt` will look something like this:


```
SLC6A4               5  1.0000
SLC6A3               4  0.8000
SLC6A2               4  0.8000
SIGMAR1              4  0.8000
HTR2C                4  0.8000
NET                  3  0.6000
KCNH2                3  0.6000
HTR2B                3  0.6000
...
```

This file is actually 50 lines long but only the top 8 are show. It shows which are the individual targets common in all of them. As we can see, [SLC6A4](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC6A4) (the serotonin transporter) is common to all 5 molecules, followed by [SLC6A3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC6A3) and [SLC6A2](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC6A2) the dopamine and noradrenaline transporter, respectively. Also the the [Sigma1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SIGMA1) Receptor and the [Serotonin 2C receptor](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR2C) are known to be affected by 4 out of 5 molecules

On the other hand we can look at `anti-dep_anno_count_T1.txt` for potential targets, as these refer to the activity profiles of similar molecules (in this case 80% similar) 

```
SLC6A4              90  0.4918
SLC6A2              38  0.2077
SLC6A3              28  0.1530
ADRA1A              16  0.0874
NET                 13  0.0710
SIGMAR1             12  0.0656
HTR2A               11  0.0601
DRD2                10  0.0546
...
```
of the 183 molecules found, the top spots were essentially the same, but the [ADRA1A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ADRA1A) Adrenergic A1 receptor and the [NET](https://www.genecards.org/cgi-bin/carddisp.pl?gene=NET), the norepinephrine transporter, appear with a high representativity. Troubling might be the presence of [HTR2A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=HTR2A), in 6% of all matches. 

To have a better look at which of the base molecules are associated, we may run the above command again with the following option `-report AD` This will create two more files that stand for the individual targets each molecule was active on. The format is a standard 'transaction type' file with one line for each molecule, followed by the found targets, named similarly to the other two files created:

* `anti-dep_annotations_T0.txt` - binding targets to each molecule in the input file
* `anti-dep_annotations_T1.txt` - binding targets to each molecule similar to the molecules in the input file.

Therefore `anti-dep_annotations_T0.txt` with the above options will look like this (end of lines arranged for this report, in the original, only one line per molecule is produced):

```
Sertraline: CYP2D6 HTR2B  CACNA1C CHRM2  SIGMAR1 CYP2C19 ADRA2C HTR2C  SLC6A3 CHRM4  ADRA2B CYP3A4
            SLC6A2 KCNH2  SLC6A4 ADRA1B KMT2A  CHRM1  HTR2A  MC5R   ADRA2A CHRM5  
Citalopran: KCNH2  SLC22A1 HTR2B  HRH1   SLC6A4 ADRA1A HTR2C  NET    SLC6A3 ADRA1B SMN2   SIGMAR1 
            ADRA1D SLC6A2 
Fluoxetine: CYP2D6 HTR3A  CACNA1C LMNA   SIGMAR1 CYP2C19 CHRM3  CYP2C9 HRH1   HTR2C  SLC6A3 ADRA2B 
            KCNK9  CYP3A4 SLC6A2 KCNH2  SLCO1B3 DRD2   SLC6A4 ACHE   CYP2B6 SLCO2B1 CYP1A2 HRH3   
            CHRM1  HTR6   HTR2A  NET    ABCB1  CYP2C8 SLCO1B1 KCNK2  ADRA2A CHRM5  
Tradozone: HTR1B  DRD2   HRH1   ADRA2C DRD3   ADRA1A HTR2B  HTR2C  HTR1A  SLC6A4 ALB    ADRA1B ADRA2A
           ADRA1D SIGMAR1 ADRA2B HTR2A  FAAH   
Venlafaxine: SLC6A4 ADRA1A SLC6A2 SLC6A3 LMNA   NET    
```

Similarly `anti-dep_annotations_T1.txt` with the above options will look like this:


```
Sertraline: CYP2D6 CYP2C9 SLC6A4 NET    SLC6A3 SLC6A2 
Citalopran: CYP2C9 SLC6A4 NFKB1  NET    SLC6A3 CYP2C19 ALOX15 CHRM1  CYP3A4 SLC6A2 
Fluoxetine: CYP2D6 HRH3   KCNH2  TSHR   SLC6A4 RORC   CACNA1C NFKB1  KCNK2  SLC6A3 KCNC1  CYP1A2 
            CYP2C19 CHRM1  CYP3A4 SLC6A2 
Tradozone: DRD2   HRH1   ADRA1A HTR1A  PDE10A SMN2   SIGMAR1 HTR2A  SLC6A2 
Venlafaxine: SLC6A4 ADRA1A SLC6A2 SLC6A3 NET
```

As an example of a possible conclusion is that in some of the 80% similar molecules to Sertraline some of them were verified active on NET and CYP2C9, strongly suggesting that Sertraline itself, although not registered as active in these two targets, might actually be so. Other types of observations might be held for other molecules.
