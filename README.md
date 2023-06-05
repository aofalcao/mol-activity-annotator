# mol-activity-annotator

`csanno` is a molecular annotator that uses ChEMBL to search to identify plausible targets. The question we are trying to solve in this program could be better defined with a simple example. Imagine we have one molecule and we want to identify probable targets in which it may be active. There are two possible situations: 

1. The molecule is present in ChEMBL and therefore all the targets for which it was found active are retrieved
2. The molecule is not found in ChEMBL. And in that case what it will do is that it will look for similar molecules (using ChEMBL's similarity search) and will output the found targets and the number of molecule similars that were actually found active there

The tool also has a input file, option where for instance we might group a set of molecules with similar physiological activity, to try to identify common or most likely targets. All of this would be actually possible to do manually using ChEMBL online, but it would be tedious and extremely prone to errors. 


`csanno` is heavily reliant on [ChEMBL web services](https://pubmed.ncbi.nlm.nih.gov/28602100/), therefore it is fundamental to have a good internet connection. Alternatively if a local ChEMBL installation is available it can be easily used with increased performance benefits.

The annotation profile of a molecule is essentially a report of the evidence that ChEMBL has for that specific structure or, if required, of the molecules that are very similar to it.
