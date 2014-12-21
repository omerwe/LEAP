import numpy as np
from optparse import OptionParser
import scipy.linalg as la
import time
import leapUtils
import scipy.linalg.blas as blas
np.set_printoptions(precision=3, linewidth=200)

parser = OptionParser()
parser.add_option('--bfilesim', metavar='bfilesim', default=None, help='Binary plink file')
parser.add_option('--extractSim', metavar='extractSim', default=None, help='SNPs subset to use')
parser.add_option('--out', metavar='out', default=None, help='output file')
parser.add_option('--pheno', metavar='pheno', default=None, help='Phenotypes file (optional), only used for identifying unphenotyped individuals')
parser.add_option('--missingPhenotype', metavar='missingPhenotype', default='-9', help='identifier for missing values (default: -9)')
(options, args) = parser.parse_args()


def eigenDecomposeMain(bfilesim, outFile, extractSim=None, pheno=None, missingPhenotype='-9'):
	#Read input files
	bed, _ = leapUtils.loadData(bfilesim, extractSim, pheno, missingPhenotype, loadSNPs=True)

	#Compute kinship matrix
	t0 = time.time()
	print 'Computing kinship matrix...'	
	XXT = leapUtils.symmetrize(blas.dsyrk(1.0, bed.val, lower=1)) / bed.val.shape[1]
	print 'Done in %0.2f'%(time.time()-t0), 'seconds'

	#Compute eigendecomposition
	S,U = leapUtils.eigenDecompose(XXT)	
	np.savez_compressed(outFile, arr_0=U, arr_1=S, XXT=XXT)

if __name__ == '__main__':
	if (options.bfilesim is None): raise Exception('bfilesim must be supplied')
	if (options.out is None):   raise Exception('output file name must be supplied')
	eigenDecomposeMain(options.bfilesim, options.out, options.extractSim, options.pheno, options.missingPhenotype)

