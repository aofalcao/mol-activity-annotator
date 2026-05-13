## Descrition of the data sets and tests


### SAR Files

[this is almost a copy from the README.md but important, as the program can read data from files only in SAR format]

SAR files are text files that should store activity data for a given target for a set of molecules and each line represents one molecule. Columns are separated by `<tab>`. The basic structur of a SAR file is

* The first column is an alphanumeric identifier for the molecule (should be unique);
* The second column contains the activity registered for that molecule (may be binary qualitative)
* The third column contains the molecule SMILES

Here is a sample SAR file, where the identifiers are the ChEMBL IDs, the middle column is the default "unknown" or irrelevant ($NA$) and the 3rd column the corresponding SMILES

```
Sertraline	NA	ClC1=CC=C([C@H]2C3=C([C@H](CC2)NC)C=CC=C3)C=C1Cl
Citalopran	NA	Fc1ccc(cc1)C3(OCc2cc(C#N)ccc23)CCCN(C)C
Fluoxetine	NA	CNCCC(c1ccccc1)Oc2ccc(cc2)C(F)(F)F
Tradozone	NA	Clc4cccc(N3CCN(CCCN1/N=C2/C=C\C=C/N2C1=O)CC3)c4
Venlafaxine	NA	OC2(C(c1ccc(OC)cc1)CN(C)C)CCCCC2
```


### Single Molecules

[examples with single molecules]

### Antidepressants

### Anticonvulsants

### Antihypertensives
