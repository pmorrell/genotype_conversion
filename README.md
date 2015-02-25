# A Python script for genotype conversion

The script takes genotype calls from [Alchemy] (http://alchemy.sourceforge.net) in the format 'AA', 'AC', 'CA', or 'CC' and converts them to a format used by various programs where the genotype is represented by the count of number of observations of the minor allele. In this code, missing data is recorded as 'NA' and is returned in that same form. The 'NA' values do not contribute to minor allele counts.

Example datasets are provided in test.txt, test2.txt, and test3.txt.


##Known Issues

###Resolved Issues
- Loci that include only 'NA' values. These genotypes are simply returned as NA. 

###Unresolved Issues
- Genotypes in some data files are represented as 'AA', 'AB', or 'BB'. This should not occur with data directly from Alchemy, but is a problem in some of our data files.
- Datafiles that in a newline '\n' at the end of the file don't work. The genotype_conversion.py returns only the headers.



