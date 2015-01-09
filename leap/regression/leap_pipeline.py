import pandas as pd
import numpy as np
import leap.leapUtils as leapUtils
import leap.leapMain as leapMain
from pysnptools.snpreader.bed import Bed

#Define analysis data
bfile = 'dataset1/dataset1'
phenoFile = bfile+'.phe'
chromosomes = xrange(1,11)
prevalence = 0.001

#Find individuals to exclude to eliminate relatedness (kinship coeff > 0.05)
bed = Bed(bfile).read().standardize()
indsToKeep = leapUtils.findRelated(bed, cutoff=0.05)

#Iterate over each chromosome
frame_list = []
for chrom in chromosomes:
	print
	print 'Analyzing chromosome', chrom, '...'

	#Create a bed object excluding SNPs from the current chromosome
	bedExclude = leapUtils.getExcludedChromosome(bfile, chrom)
	
	#Create a bed object including only SNPs from the current chromosome
	bedTest = leapUtils.getChromosome(bfile, chrom)	
	
	#Compute eigendecomposition for the data
	eigenFile = 'temp_eigen.npz'
	eigen = leapMain.eigenDecompose(bedExclude, outFile=eigenFile)
	
	#compute heritability explained by this data
	h2 = leapMain.calcH2(phenoFile, prevalence, eigen, keepArr=indsToKeep, h2coeff=1.0)
	
	#Compute liabilities explained by this data
	liabs = leapMain.probit(bedExclude, phenoFile, h2, prevalence, eigen, keepArr=indsToKeep)
	
	#perform GWAS, using the liabilities as the observed phenotypes
	results_df = leapMain.leapGwas(bedExclude, bedTest, liabs, h2)
	frame_list.append(results_df)

#Join together the results of all chromosomes, and print the top ranking SNPs
frame = pd.concat(frame_list)
frame.sort("PValue", inplace=True)
frame.index = np.arange(len(frame))
print 'Top 10 most associated SNPs:'
print frame.head(n=10)


	
	
	