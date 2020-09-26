import numpy as np
import argparse
import time
import scipy.linalg.blas as blas
from . import leapUtils
np.set_printoptions(precision=3, linewidth=200)
from . import leapMain


def findRelated(bed, outFile, cutoff, kinshipFile=None):

	bed = leapUtils._fixupBed(bed)
	
	keepArr = leapUtils.findRelated(bed, cutoff, kinshipFile)
	if (outFile is not None):
		print('Printing output to', outFile)
		f = open(outFile, 'w')
		for i, (fid,iid) in enumerate(bed.iid):
			if (keepArr[i]): f.write(fid + ' ' + iid + ' 0\n')
			else: f.write(fid + ' ' + iid + ' 1\n')
		f.close()
	return keepArr
	
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('--bfilesim', metavar='bfilesim', default=None, help='Binary plink file')
	parser.add_argument('--kinship', metavar='kinsip', default=None, help='A kinship matrix represented in a text file. Note that this matrix must correspond exactly to the phenotypes file, unlike the bfilesim file option.')
	parser.add_argument('--extractSim', metavar='extractSim', default=None, help='extractSim file')
	parser.add_argument('--cutoff', metavar='cutoff', type=float, default=0.05, help='Relationship cutoff (default 0.05)')
	parser.add_argument('--out', metavar='out', default=None, help='output file')

	parser.add_argument('--pheno', metavar='pheno', default=None, help='Phenotypes file (optional), only used for identifying unphenotyped individuals')
	parser.add_argument('--missingPhenotype', metavar='missingPhenotype', default='-9', help='identifier for missing values (default: -9)')
	args = parser.parse_args()

	if (args.bfilesim is None and args.kinship is None): raise Exception('bfilesim or kinship must be supplied')
	if (args.bfilesim is not None and args.kinship is not None): raise Exception('bfilesim and kinship cannot both be supplied')
	if (args.out is None):   raise Exception('output file name must be supplied')
	if (args.bfilesim is not None): bed, _ = leapUtils.loadData(args.bfilesim, args.extractSim, args.pheno, args.missingPhenotype, loadSNPs=True, standardize=True)
	else: bed=None
	leapMain.findRelated(bed, args.out, args.cutoff, args.kinship)
	
