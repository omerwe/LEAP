import numpy as np
import argparse
import scipy.linalg as la
import time
import leapUtils
import scipy.linalg.blas as blas
import leapMain
np.set_printoptions(precision=3, linewidth=200)

def eigenDecompose(bed, outFile=None):

	bed = leapUtils._fixupBed(bed)

	#Compute kinship matrix
	t0 = time.time()
	print 'Computing kinship matrix...'	
	XXT = leapUtils.symmetrize(blas.dsyrk(1.0, bed.val, lower=1)) / bed.val.shape[1]
	print 'Done in %0.2f'%(time.time()-t0), 'seconds'

	#Compute eigendecomposition
	S,U = leapUtils.eigenDecompose(XXT)
	if (outFile is not None): np.savez_compressed(outFile, arr_0=U, arr_1=S, XXT=XXT)	
	eigen = dict([])
	eigen['XXT'] = XXT
	eigen['arr_0'] = U
	eigen['arr_1'] = S
	return eigen

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--bfilesim', metavar='bfilesim', default=None, help='Binary plink file')
	parser.add_argument('--extractSim', metavar='extractSim', default=None, help='SNPs subset to use')
	parser.add_argument('--out', metavar='out', default=None, help='output file')
	parser.add_argument('--pheno', metavar='pheno', default=None, help='Phenotypes file (optional), only used for identifying unphenotyped individuals')
	parser.add_argument('--missingPhenotype', metavar='missingPhenotype', default='-9', help='identifier for missing values (default: -9)')
	args = parser.parse_args()
	
	if (args.bfilesim is None): raise Exception('bfilesim must be supplied')
	if (args.out is None):   raise Exception('output file name must be supplied')
	
	#Read input files
	bed, _ = leapUtils.loadData(args.bfilesim, args.extractSim, args.pheno, args.missingPhenotype, loadSNPs=True)
	leapMain.eigenDecompose(bed, args.out)

