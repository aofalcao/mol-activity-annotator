## Descrition of the data sets and tests


### SAR Files

[this is almost a copy from the README.md but important, as the program can read data from files only in SAR format]

SAR files are text files that should store activity data for a given target for a set of molecules and each line represents one molecule. Columns are separated by `<tab>`. The basic structur of a SAR file is

* The first column is an alphanumeric identifier for the molecule (should be unique);
* The second column contains the activity registered for that molecule (may be binary qualitative)
* The third column contains the molecule SMILES

Here is a sample SAR file, where the identifiers are the ChEMBL IDs, the middle column is the default "unknown" or irrelevant ($NA$) and the 3rd column the corresponding SMILES

```
CHEMBL158973	NA	CN(C)CCSC(C)(C)C
CHEMBL476516	NA	NCCc1ccc(Br)cc1
CHEMBL309689	NA	Oc1noc2c1CCNCC2
```


### Single Molecules

[examples with single molecules]

### Antidepressants

### Anticonvulsants

### Antihypertensives
