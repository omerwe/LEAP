import numpy as np
from optparse import OptionParser
import scipy.linalg as la
import scipy.stats as stats
import scipy.linalg.blas as blas
import csv
import sklearn.linear_model
import time
import sys
import VertexCut as vc
from pysnptools.pysnptools.snpreader.bed import Bed
from fastlmm.pyplink import plink
import pysnptools.pysnptools.util.util as pyutil
np.set_printoptions(precision=3, linewidth=200)
import leapUtils

parser = OptionParser()
parser.add_option('--bfilesim', metavar='bfilesim', default=None, help='Binary plink file')
parser.add_option('--extractSim', metavar='extractSim', default=None, help='SNPs subset to use')
parser.add_option('--prev', metavar='prev', type=float, default=None, help='Trait prevalence')
parser.add_option('--numRemovePCs', metavar='numRemovePCs', type=int, default=10, help='Number of principal components to fit')
parser.add_option('--pheno', metavar='pheno', default=None, help='Phenotype file in Plink format')
parser.add_option('--eigen', metavar='eigen', default=None, help='eigen file')
parser.add_option('--related', metavar='related', default=None, help='relatedness file')

parser.add_option('--lowtail', metavar='lowtail', type=int, default=0, help='Assume that both tails of the liabilities distribution are oversampled (0 or 1 - default 0)')
parser.add_option('--h2coeff', metavar='h2coeff', type=float, default=0.85, help='Heritability coefficient (set to 1.0 for synthetic data)')
parser.add_option('--relCutoff', metavar='relCutoff', type=float, default=0.05, help='relatedness cutoff (set to negative value to override relatedness check)')
parser.add_option('--missingPhenotype', metavar='missingPhenotype', type=float, default=-9, help='identifier for missing values (default: -9)')
(options, args) = parser.parse_args()

if (options.bfilesim is None): raise Exception('--bfilesim must be supplied')
if (options.prev is None): raise Exception('--prev must be supplied')
if (options.pheno is None): raise Exception('--pheno must be supplied')


def calcLiabThreholds(U, s, keepArr, phe):

	#Run logistic regression
	G = U[:, :options.numRemovePCs] * np.sqrt(s[:options.numRemovePCs])
	Logreg = sklearn.linear_model.LogisticRegression(penalty='l2', C=500000, fit_intercept=True)
	Logreg.fit(G[keepArr, :options.numRemovePCs], phe[keepArr])

	#Compute individual thresholds
	Pi = Logreg.predict_proba(G[:, :options.numRemovePCs])[:,1]

	#Compute thresholds and save to files
	P = np.sum(phe==1) / float(phe.shape[0])
	K = options.prev
	Ki = K*(1-P) / (P*(1-K)) * Pi / (1 + K*(1-P) / (P*(1-K))*Pi - Pi)
	thresholds = stats.norm(0,1).isf(Ki)
	thresholds[Ki>=1.] = -999999999
	thresholds[Ki<=0.] = 999999999
	
	return Pi, thresholds
	
	
def calcH2Continuous_twotails(XXT, phe, keepArr):

	print 'computing h2 for a two-tails ascertained study...'
	
	XXT = XXT[np.ix_(keepArr, keepArr)]
	phe = phe[keepArr]	

	t1 = stats.norm(0,1).ppf(options.prev)
	t2 = stats.norm(0,1).isf(options.prev)
	phit1 = stats.norm(0,1).pdf(t1)
	phit2 = stats.norm(0,1).pdf(t2)
	
	K1 = options.prev
	K2 = options.prev
	
	xCoeff = ((phit2*t2 - phit1*t1 + K1 + K2)**2 * (K1+K2)**2 - (phit2-phit1)**4) / (K1 + K2)**4
	intersect = ((phit2-phit1) / (K1+K2))**2
		
	pheMean = 0
	pheVar = 1
	
	x = (xCoeff * options.h2coeff) * XXT 
	y = np.outer((phe-pheMean)/np.sqrt(pheVar), (phe-pheMean)/np.sqrt(pheVar))
	y -= intersect
	
	y = y[np.triu_indices(y.shape[0], 1)]
	x = x[np.triu_indices(x.shape[0], 1)]
	
	slope, intercept, rValue, pValue, stdErr = stats.linregress(x,y)	
	return slope
	
	

def calcH2Continuous(XXT, phe, keepArr):
	t = stats.norm(0,1).isf(options.prev)
	phit = stats.norm(0,1).pdf(t)
	
	K1 = 1 - options.prev
	K2 = 1 - K1
	P = np.sum(phe<t) / float(phe.shape[0])	
	P2 = 1.0
	P1 = K2*P2*P / (K1*(1-P))
	R = P2 / P1
	
	XXT = XXT[np.ix_(keepArr, keepArr)]
	phe = phe[keepArr]
	
	xCoeff = (((R-1)*phit*t + K1 + R*K2)**2 * (K1+R*K2)**2 - ((R-1)*phit)**4) / (K1 + R*K2)**4
	x = (xCoeff * options.h2coeff) * XXT 
	pheMean = 0
	pheVar = 1	
	y = np.outer((phe-pheMean) / np.sqrt(pheVar), (phe-pheMean)/np.sqrt(pheVar))
	y -= ((R-1)*phit / (K1+R*K2))**2
	
	y = y[np.triu_indices(y.shape[0], 1)]
	x = x[np.triu_indices(x.shape[0], 1)]
	
	slope, intercept, rValue, pValue, stdErr = stats.linregress(x,y)
	return slope
	
	
	
def calcH2Binary(XXT, phe, probs, thresholds, keepArr):
	K = options.prev
	P = np.sum(phe>0) / float(phe.shape[0])
	
	XXT = XXT[np.ix_(keepArr, keepArr)]
	phe = phe[keepArr]
	
	if (thresholds is None):
		t = stats.norm(0,1).isf(K)
		phit = stats.norm(0,1).pdf(t)
		xCoeff = P*(1-P) / (K**2 * (1-K)**2) * phit**2 * options.h2coeff
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
		x = (Atag0 / B0 * options.h2coeff) * XXT	
	
	y = y[np.triu_indices(y.shape[0], 1)]
	x = x[np.triu_indices(x.shape[0], 1)]
	
	slope, intercept, rValue, pValue, stdErr = stats.linregress(x,y)
	return slope
		
######## main ######

#Determine if this is a case-control study
pheOrig = np.loadtxt(options.pheno, usecols=[2])
pheOrig = pheOrig[pheOrig != options.missingPhenotype]
pheUnique = np.unique(pheOrig)
if (pheUnique.shape[0] < 2): raise Exception('Less than two different phenotypes observed')
isCaseControl = (pheUnique.shape[0] == 2)

#Read eigen file	
if (options.eigen is not None): eigen = np.load(options.eigen)

#Read bfilesim and pheno file for heritability computation
loadSNPs = (options.eigen is None)	
bed, phe = leapUtils.loadData(options.bfilesim, options.extractSim, options.pheno, options.missingPhenotype, loadSNPs=loadSNPs, standardize=True)

# # # #Standardize SNPs according to frq file
# # # pheThreshold = (np.mean(pheUnique) if isCaseControl else stats.norm(0,1).isf(options.prev))
# # # empMean = np.mean(bed.val, axis=0) / 2.0
# # # bed.val[:, empMean>0.5] = 2 - bed.val[:, empMean>0.5]
# # # controls = (phe<pheThreshold)
# # # mafs = bed.val[controls, :].mean(axis=0)/2.0
# # # #mafs = np.loadtxt(options.bfilesim+'.frq', usecols=[1,2]).mean(axis=1)
# # # bed.val -= 2*mafs
# # # bed.val /= np.sqrt(2*mafs*(1-mafs))


#Compute kinship matrix
if (options.eigen is None): XXT = leapUtils.symmetrize(blas.dsyrk(1.0, bed.val, lower=1)) / bed.val.shape[1]
else: XXT = eigen['XXT'] / bed.sid.shape[0]

#Compute relatedness
if (options.relCutoff <= 0): keepArr = np.ones(bed.iid.shape[0], dtype=bool)
else:
	bed2 = bed
	if (options.related is None):
		if (options.extractSim is not None): bed2, _ = leapUtils.loadData(options.bfilesim, None, options.pheno, options.missingPhenotype, loadSNPs=True)			
		keepArr = leapUtils.findRelated(bed2, options.relCutoff)
	else:
		keepArr = leapUtils.loadRelatedFile(bed, options.related)	

#Remove top PCs from kinship matrix
if (options.numRemovePCs > 0):
	if (options.eigen is None): s,U = leapUtils.eigenDecompose(XXT)
	else: s, U = eigen['s']**2 / bed.sid.shape[0], eigen['U']	
	
	XXT -= (U[:, :options.numRemovePCs]*s[:options.numRemovePCs]).dot(U[:, :options.numRemovePCs].T)
	
if isCaseControl:
	print 'Computing h2 for a binary phenotype'
	pheMean = phe.mean()	
	phe[phe <= pheMean] = 0
	phe[phe > pheMean] = 1	
	if (options.numRemovePCs > 0):
		probs, thresholds = calcLiabThreholds(U, s, keepArr, phe)
		h2 = calcH2Binary(XXT, phe, probs, thresholds, keepArr)
	else: h2 = calcH2Binary(XXT, phe, None, None, keepArr)
else:
	print 'Computing h2 for a continuous phenotype'
	if (options.lowtail == 0): h2 = calcH2Continuous(XXT, phe, keepArr)
	else: h2 = calcH2Continuous_twotails(XXT, phe, keepArr)
	
if (h2 <= 0): raise Exception("Negative heritability found. Exitting...")	
if (np.isnan(h2)): raise Exception("Invalid heritability estimate. Please double-check your input for any errors.")	
	
print 'h2: %0.6f'%h2





