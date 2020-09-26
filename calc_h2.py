import numpy as np
import argparse
import scipy.stats as stats
import scipy.linalg.blas as blas
import sklearn.linear_model
import time
import sys
np.set_printoptions(precision=3, linewidth=200)
from . import leapUtils
from . import leapMain


def calcLiabThreholds(U, S, keepArr, phe, numRemovePCs, prev, covar):

	#Run logistic regression
	if (numRemovePCs > 0):
		G = U[:, -numRemovePCs:] * np.sqrt(S[-numRemovePCs:])
	else:
		G = np.empty((phe.shape[0], 0))
	if (covar is not None): G = np.concatenate((G, covar), axis=1)
	
	Logreg = sklearn.linear_model.LogisticRegression(penalty='l2', C=500000, fit_intercept=True)
	Logreg.fit(G[keepArr, :], phe[keepArr])
	
	#Compute individual thresholds
	Pi = Logreg.predict_proba(G)[:,1]

	#Compute thresholds
	P = np.sum(phe==1) / float(phe.shape[0])
	K = prev
	Ki = K*(1-P) / (P*(1-K)) * Pi / (1 + K*(1-P) / (P*(1-K))*Pi - Pi)
	thresholds = stats.norm(0,1).isf(Ki)
	thresholds[Ki>=1.] = -999999999
	thresholds[Ki<=0.] = 999999999
	
	return Pi, thresholds
	
	
def calcH2Continuous_twotails(XXT, phe, keepArr, prev):

	print('computing h2 for a two-tails ascertained study...')
	
	XXT = XXT[np.ix_(keepArr, keepArr)]
	phe = phe[keepArr]	

	t1 = stats.norm(0,1).ppf(prev)
	t2 = stats.norm(0,1).isf(prev)
	phit1 = stats.norm(0,1).pdf(t1)
	phit2 = stats.norm(0,1).pdf(t2)
	
	K1 = prev
	K2 = prev
	
	xCoeff = ((phit2*t2 - phit1*t1 + K1 + K2)**2 * (K1+K2)**2 - (phit2-phit1)**4) / (K1 + K2)**4
	intersect = ((phit2-phit1) / (K1+K2))**2
		
	pheMean = 0
	pheVar = 1
	
	x = xCoeff * XXT 
	y = np.outer((phe-pheMean)/np.sqrt(pheVar), (phe-pheMean)/np.sqrt(pheVar))
	y -= intersect
	
	y = y[np.triu_indices(y.shape[0], 1)]
	x = x[np.triu_indices(x.shape[0], 1)]
	
	slope, intercept, rValue, pValue, stdErr = stats.linregress(x,y)	
	return slope
	
	

def calcH2Continuous(XXT, phe, keepArr, prev):
	t = stats.norm(0,1).isf(prev)
	phit = stats.norm(0,1).pdf(t)
	
	K1 = 1 - prev
	K2 = 1 - K1
	P = np.sum(phe<t) / float(phe.shape[0])	
	P2 = 1.0
	P1 = K2*P2*P / (K1*(1-P))
	R = P2 / P1
	
	XXT = XXT[np.ix_(keepArr, keepArr)]
	phe = phe[keepArr]
	
	xCoeff = (((R-1)*phit*t + K1 + R*K2)**2 * (K1+R*K2)**2 - ((R-1)*phit)**4) / (K1 + R*K2)**4
	x = xCoeff * XXT 
	pheMean = 0
	pheVar = 1	
	y = np.outer((phe-pheMean) / np.sqrt(pheVar), (phe-pheMean)/np.sqrt(pheVar))
	y -= ((R-1)*phit / (K1+R*K2))**2
	
	y = y[np.triu_indices(y.shape[0], 1)]
	x = x[np.triu_indices(x.shape[0], 1)]
	
	slope, intercept, rValue, pValue, stdErr = stats.linregress(x,y)
	return slope
	
	
	
def calcH2Binary(XXT, phe, probs, thresholds, keepArr, prev):
	K = prev
	P = np.sum(phe>0) / float(phe.shape[0])
	
	XXT = XXT[np.ix_(keepArr, keepArr)]
	phe = phe[keepArr]
	
	if (thresholds is None):
		t = stats.norm(0,1).isf(K)
		phit = stats.norm(0,1).pdf(t)
		xCoeff = P*(1-P) / (K**2 * (1-K)**2) * phit**2
		y = np.outer((phe-P) / np.sqrt(P*(1-P)), (phe-P) / np.sqrt(P*(1-P)))
		x = xCoeff * XXT
		
	else:
		probs = probs[keepArr]
		thresholds = thresholds[keepArr]
		Ki = K*(1-P) / (P*(1-K)) * probs / (1 + K*(1-P) / (P*(1-K))*probs - probs)
		phit = stats.norm(0,1).pdf(thresholds)	
		probsInvOuter = np.outer(probs*(1-probs), probs*(1-probs))
		y = np.outer(phe-probs, phe-probs) / np.sqrt(probsInvOuter)	
		sumProbs = np.tile(np.column_stack(probs).T, (1,probs.shape[0])) + np.tile(probs, (probs.shape[0], 1))
		Atag0 = np.outer(phit, phit) * (1 - (sumProbs)*(P-K)/(P*(1-K)) + np.outer(probs, probs)*(((P-K)/(P*(1-K)))**2)) / np.sqrt(probsInvOuter)
		B0 = np.outer(Ki + (1-Ki)*(K*(1-P))/(P*(1-K)), Ki + (1-Ki)*(K*(1-P))/(P*(1-K)))
		x = Atag0 / B0 * XXT	
	
	y = y[np.triu_indices(y.shape[0], 1)]
	x = x[np.triu_indices(x.shape[0], 1)]
	
	slope, intercept, rValue, pValue, stdErr = stats.linregress(x,y)
	return slope
		
		
		
def calc_h2(pheno, prev, eigen, keepArr, covar, numRemovePCs, lowtail):

	pheno = leapUtils._fixup_pheno(pheno)

	#Extract phenotype
	if isinstance(pheno, dict):	phe = pheno['vals']	
	else: phe = pheno		
	if (len(phe.shape)==2):
		if (phe.shape[1]==1): phe=phe[:,0]
		else: raise Exception('More than one phenotype found')		
	if (keepArr is None): keepArr = np.ones(phe.shape[0], dtype=np.bool)
	
	#Compute kinship matrix	
	XXT = eigen['XXT']
	
	#Remove top PCs from kinship matrix
	if (numRemovePCs > 0):
		if (eigen is None): S,U = leapUtils.eigenDecompose(XXT)
		else: S, U = eigen['arr_1'], eigen['arr_0']		
		print('Removing the top', numRemovePCs, 'PCs from the kinship matrix')
		XXT -= (U[:, -numRemovePCs:]*S[-numRemovePCs:]).dot(U[:, -numRemovePCs:].T)		
	else:
		U, S = None, None
		
	#Determine if this is a case-control study
	pheUnique = np.unique(phe)
	if (pheUnique.shape[0] < 2): raise Exception('Less than two different phenotypes observed')
	isCaseControl = (pheUnique.shape[0] == 2)
	
	if isCaseControl:
		print('Computing h2 for a binary phenotype')
		pheMean = phe.mean()	
		phe[phe <= pheMean] = 0
		phe[phe > pheMean] = 1
		if (numRemovePCs > 0 or covar is not None):
			probs, thresholds = calcLiabThreholds(U, S, keepArr, phe, numRemovePCs, prev, covar)
			h2 = calcH2Binary(XXT, phe, probs, thresholds, keepArr, prev)
		else: h2 = calcH2Binary(XXT, phe, None, None, keepArr, prev)
	else:
		if (covar is not None): raise Exception('Covariates with a continuous phenotype are currently not supported')
		print('Computing h2 for a continuous phenotype')
		if (not lowtail): h2 = calcH2Continuous(XXT, phe, keepArr, prev)
		else: h2 = calcH2Continuous_twotails(XXT, phe, keepArr, prev)
		
	if (h2 <= 0): raise Exception("Negative heritability found. Exitting...")	
	if (np.isnan(h2)): raise Exception("Invalid heritability estimate. Please double-check your input for any errors.")	
		
	print('h2: %0.6f'%h2)
	return h2


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('--bfilesim', metavar='bfilesim', default=None, help='Binary plink file')
	parser.add_argument('--extractSim', metavar='extractSim', default=None, help='SNPs subset to use')
	parser.add_argument('--prev', metavar='prev', type=float, default=None, help='Trait prevalence')
	parser.add_argument('--numRemovePCs', metavar='numRemovePCs', type=int, default=10, help='Number of principal components to fit')
	parser.add_argument('--pheno', metavar='pheno', default=None, help='Phenotype file in Plink format')
	parser.add_argument('--eigen', metavar='eigen', default=None, help='eigen file')
	parser.add_argument('--related', metavar='related', default=None, help='relatedness file')
	parser.add_argument('--covar', metavar='covar', default=None, help='covariates file')

	parser.add_argument('--lowtail', metavar='lowtail', type=int, default=0, help='Assume that both tails of the liabilities distribution are oversampled (0 or 1 - default 0)')
	parser.add_argument('--relCutoff', metavar='relCutoff', type=float, default=0.05, help='relatedness cutoff (set to negative value to override relatedness check)')
	parser.add_argument('--missingPhenotype', metavar='missingPhenotype', default='-9', help='identifier for missing values (default: -9)')
	args = parser.parse_args()


	if (args.bfilesim is None): raise Exception('--bfilesim must be supplied')
	if (args.prev is None): raise Exception('--prev must be supplied')
	if (args.pheno is None): raise Exception('--pheno must be supplied')
	
	#Read bfilesim and pheno file for heritability computation	
	bed, phe = leapUtils.loadData(args.bfilesim, args.extractSim, args.pheno, args.missingPhenotype, loadSNPs=(args.eigen is None), standardize=True)
	
	#Read/create eigendecomposition
	if (args.eigen is not None): eigen = np.load(args.eigen)
	else:
		from . import eigenDecompose
		eigen = eigenDecompose.eigenDecompose(bed)	

	#Compute relatedness
	if (args.relCutoff <= 0): keepArr = np.ones(bed.iid.shape[0], dtype=bool)
	else:		
		if (args.related is None):
			bed2 = bed
			if (args.extractSim is not None or args.eigen is not None): bed2, _ = leapUtils.loadData(args.bfilesim, None, args.pheno, args.missingPhenotype, loadSNPs=True)			
			keepArr = leapUtils.findRelated(bed2, args.relCutoff)
		else:
			keepArr = leapUtils.loadRelatedFile(bed, args.related)
			
			
	#Read covar file
	if (args.covar is not None):		
		covar = leapUtils.loadCovars(bed, args.covar)	
		covar -= covar.mean()
		covar /= covar.std()
		print('Read', covar.shape[1], 'covariates from file')
	else:
		covar = None		
	
	leapMain.calcH2(phe, args.prev, eigen, keepArr, covar, args.numRemovePCs, args.lowtail==1)





