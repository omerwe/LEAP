import numpy as np
import time
import sys
import argparse
import leapUtils
import leapMain
import fastlmm.association
np.set_printoptions(precision=4, linewidth=200)


def gwas(bedSim, bedTest, pheno, h2, outFile, eigenFile, covar):

	bedSim, pheno = leapUtils._fixupBedAndPheno(bedSim, pheno)
	bedTest, pheno = leapUtils._fixupBedAndPheno(bedTest, pheno)

	#Run GWAS	
	logdelta = np.log(1.0/h2 - 1)
	G0 = (bedSim if eigenFile is None else None)
	print 'Performing LEAP GWAS...'
	results_df = fastlmm.association.single_snp(bedTest, pheno, G0=G0, covar=covar, output_file_name=outFile, log_delta=logdelta, cache_file=eigenFile)
	return results_df
	

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('--bfilesim', metavar='bfilesim', default=None, help='Binary plink file')
	parser.add_argument('--bfile', metavar='bfile', default=None, help='Binary plink file to test')
	parser.add_argument('--pheno', metavar='pheno', default=None, help='Phenotype file in Plink format')
	parser.add_argument('--eigen', metavar='eigen', default=None, help='Eigendecompositon file')
	parser.add_argument('--h2', metavar='h2', type=float, default=None, help='h2 value')
	parser.add_argument('--extractSim', metavar='extractSim', default=None, help='SNPs subset to use')
	parser.add_argument('--extract', metavar='extract', default=None, help='SNPs subset to test')
	parser.add_argument('--out', metavar='out', default=None, help='output file')
	parser.add_argument('--covar', metavar='covar', default=None, help='Covariates file')
	parser.add_argument('--missingPhenotype', metavar='missingPhenotype', default='-9', help='identifier for missing values (default: -9)')
	args = parser.parse_args()

	if (args.bfile is None): raise Exception('bfile must be supplied')
	if (args.bfilesim is None): raise Exception('bfilesim must be supplied')
	if (args.pheno is None): raise Exception('phenotype file must be supplied')
	if (args.out is None):   raise Exception('output file name must be supplied')
	if (args.h2 is None): raise Exception('h2 must be supplied')
	
	#Read bfile and pheno file
	bedSim, _ = leapUtils.loadData(args.bfilesim, args.extractSim, args.pheno, args.missingPhenotype, loadSNPs=True)
	bedTest, _ = leapUtils.loadData(args.bfile, args.extract, args.pheno, args.missingPhenotype, loadSNPs=True)

	#Read covariates
	if (args.covar is not None):		
		covar = leapUtils.loadCovars(bed, covar)			
		print 'Read', covarsMat.shape[1], 'covariates from file'
	else: covar = None
	
	leapMain.leapGwas(bedSim, bedTest, args.pheno, args.h2, args.out, eigenFile=args.eigen, covar=covar)

	
