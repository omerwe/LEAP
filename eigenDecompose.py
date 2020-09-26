import numpy as np
import argparse
import scipy.linalg as la
import time
from . import leapUtils
import scipy.linalg.blas as blas
from . import leapMain
np.set_printoptions(precision=3, linewidth=200)

def eigenDecompose(bed, kinshipFile=None, outFile=None, ignore_neig=False):

	if (kinshipFile is None):
		#Compute kinship matrix
		bed = leapUtils._fixupBed(bed)
		t0 = time.time()
		print('Computing kinship matrix...')	
		XXT = leapUtils.symmetrize(blas.dsyrk(1.0, bed.val, lower=1)) / bed.val.shape[1]
		print('Done in %0.2f'%(time.time()-t0), 'seconds')
	else:
		XXT = np.loadtxt(kinshipFile)

	#Compute eigendecomposition
	S,U = leapUtils.eigenDecompose(XXT, ignore_neig)
	if (outFile is not None): np.savez_compressed(outFile, arr_0=U, arr_1=S, XXT=XXT)	
	eigen = dict([])
	eigen['XXT'] = XXT
	eigen['arr_0'] = U
	eigen['arr_1'] = S
	return eigen

	
	
if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument('--bfilesim', metavar='bfilesim', default=None, help='Binary plink file')
	parser.add_argument('--kinship', metavar='kinship', default=None, help='A kinship matrix represented in a text file. Note that this matrix must correspond exactly to the phenotypes file, unlike the bfilesim file option.')
	parser.add_argument('--extractSim', metavar='extractSim', default=None, help='SNPs subset to use')
	parser.add_argument('--out', metavar='out', default=None, help='output file')
	parser.add_argument('--pheno', metavar='pheno', default=None, help='Phenotypes file (optional), only used for identifying unphenotyped individuals')
	parser.add_argument('--missingPhenotype', metavar='missingPhenotype', default='-9', help='identifier for missing values (default: -9)')
	parser.add_argument('--ignore_neig', metavar='ignore_neig', type=int, default=0, help='if set to 1, negative eigenvalues will be set to 0 and consequently ignored.')
	args = parser.parse_args()
	
	if (args.bfilesim is None and args.kinship is None): raise Exception('bfilesim or kinship must be supplied')
	if (args.bfilesim is not None and args.kinship is not None): raise Exception('bfilesim and kinship cannot both be supplied')
	if (args.out is None):   raise Exception('output file name must be supplied')
	
	#Read input files
	if (args.bfilesim is not None): bed, _ = leapUtils.loadData(args.bfilesim, args.extractSim, args.pheno, args.missingPhenotype, loadSNPs=True)
	else: bed=None
	
	leapMain.eigenDecompose(bed, kinshipFile=args.kinship, outFile=args.out, ignore_neig=args.ignore_neig>0)

