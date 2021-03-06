# A Python script for genotype conversion

The script takes genotype calls from [Alchemy] (http://alchemy.sourceforge.net) in the format 'AA', 'AC', 'CA', 'CC', CG, etc. and converts them to a format where the genotype is represented by the count of number of observations of the minor allele. Genotypes are reported as '0', '1', or '2' based on the number of minor alleles. This format is used by various programs, including the interlocus linkage disequilibrium R tool [LDcorsv] (http://cran.r-project.org/web/packages/LDcorSV/index.html).  Missing data is recorded as 'NA' and is returned in that same form. The 'NA' values do not contribute to minor allele counts.

Example datasets are provided in test.txt, test2.txt, and test3.txt.


##Known Issues

###Resolved Issues
- Loci that include only 'NA' values. These genotypes are simply returned as NA. 

###Unresolved Issues
- Genotypes in some data files are represented as 'AA', 'AB', or 'BB'. This should generally not occur with data directly from Alchemy, but is a problem in some of our data files.
- Datafiles with a newline ('\n') at the end of the file don't work. The current version of genotype_conversion.py returns only the headers.

- Potential additions to the script
	- the option to process data sets with 'AA', 'AB', 'BB' genotype calls
	- better exception handling for calls that don't match those in the defined set
	- exception handling for loci that include more than two variants at a locus
	- options to write to the '.ped' format used for [Plink] (https://www.cog-genomics.org/plink2) genotype data
        - this would also require a '.map' which includes SNP name and physical (or genetic) map position, which aren't in the Alchemy output currently processed 

#Quality control (QC) on SNP genotyping matrix

##An R script for processing a SNP genotyping matrix in minor allele count format

For allele count format, see the description above.

It is often important to perform QC on a SNP data set before proceeding with analysis. 

This code permits reading in a data set, specified at the top of the code and then filters on the criteria for heterozygosity, missingness, etc.


 
