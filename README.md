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

