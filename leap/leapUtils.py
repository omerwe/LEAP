import numpy as np
from optparse import OptionParser
import scipy.linalg as la
import scipy.stats as stats
import scipy.linalg.blas as blas
import pandas as pd
import csv
import time
import fastlmm.util.VertexCut as vc
from pysnptools.snpreader.bed import Bed
import pysnptools.util as pstutil
import pysnptools.util.pheno as phenoUtils
np.set_printoptions(precision=3, linewidth=200)



def loadData(bfile, extractSim, phenoFile, missingPhenotype='-9', loadSNPs=False, standardize=True):
	bed = Bed(bfile)
	
	if (extractSim is not None):
		f = open(extractSim)
		csvReader = csv.reader(f)
		extractSnpsSet = set([])
		for l in csvReader: extractSnpsSet.add(l[0])			
		f.close()		
		keepSnpsInds = [i for i in xrange(bed.sid.shape[0]) if bed.sid[i] in extractSnpsSet]		
		bed = bed[:, keepSnpsInds]
		
	phe = None
	if (phenoFile is not None):	bed, phe = loadPheno(bed, phenoFile, missingPhenotype)
	
	if (loadSNPs):
		bed = bed.read()
		if (standardize): bed = bed.standardize()	
	
	return bed, phe
	
	
def loadPheno(bed, phenoFile, missingPhenotype='-9', keepDict=False):
	pheno = phenoUtils.loadOnePhen(phenoFile, missing=missingPhenotype, vectorize=True)
	checkIntersection(bed, pheno, 'phenotypes')
	bed, pheno = pstutil.intersect_apply([bed, pheno])
	if (not keepDict): pheno = pheno['vals']
	return bed, pheno
	
	
def checkIntersection(bed, fileDict, fileStr, checkSuperSet=False):
	bedSet = set((b[0], b[1]) for b in bed.iid)
	fileSet = set((b[0], b[1]) for b in fileDict['iid'])
	
	if checkSuperSet:
		if (not fileSet.issuperset(bedSet)): raise Exception(fileStr + " file does not include all individuals in the bfile")
	
	intersectSet = bedSet.intersection(fileSet)
	if (len(intersectSet) != len (bedSet)):
		print len(intersectSet), 'individuals appear in both the plink file and the', fileStr, 'file'

	
def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())
	
	

def loadRelatedFile(bed, relFile):
	relatedDict = phenoUtils.loadOnePhen(relFile, vectorize=True)
	checkIntersection(bed, relatedDict, 'relatedness', checkSuperSet=True)
	_, relatedDict = pstutil.intersect_apply([bed, relatedDict])
	related = relatedDict['vals']
	keepArr = (related < 0.5)
	print np.sum(~keepArr), 'individuals will be removed due to high relatedness'
	return keepArr
	
	
def findRelated(bed, cutoff):
	print 'Computing kinship matrix...'
	t0 = time.time()	
	XXT = symmetrize(blas.dsyrk(1.0, bed.val, lower=1) / bed.val.shape[1])
	print 'Done in %0.2f'%(time.time()-t0), 'seconds'

	#Find related individuals
	removeSet = set(np.sort(vc.VertexCut().work(XXT, cutoff))) #These are the indexes of the IIDs to remove		
	print 'Marking', len(removeSet), 'individuals to be removed due to high relatedness'
	
	#keepArr = np.array([(1 if iid in keepSet else 0) for iid in bed.iid], dtype=bool)	
	keepArr = np.ones(bed.iid.shape[0], dtype=bool)
	for i in removeSet: keepArr[i] = False	
	return keepArr
	
	
	
def eigenDecompose(XXT):
	t0 = time.time()
	print 'Computing eigendecomposition...'
	s,U = la.eigh(XXT)
	if (np.min(s) < -1e-4): raise Exception('Negative eigenvalues found')
	s[s<0]=0	
	ind = np.argsort(s)
	ind = ind[s>1e-12]
	U = U[:, ind]
	s = s[ind]
	print 'Done in %0.2f'%(time.time()-t0), 'seconds'
	return s,U
	
	

def loadCovars(bed, covarFile):
	covarsDict = phenoUtils.loadPhen(covarFile)
	checkIntersection(bed, covarsDict, 'covariates', checkSuperSet=True)
	_, covarsDict = pstutil.intersect_apply([bed, covarsDict])
	covar = covarsDict['vals']
	return covar	
	
def getSNPCovarsMatrix(bed, resfile, pthresh, mindist):
	snpNameToNumDict = dict([])
	for i,s in enumerate(bed.sid): snpNameToNumDict[s] = i	

	f = open(resfile)
	csvReader = csv.reader(f, delimiter="\t")
	csvReader.next()	
	significantSNPs = []
	significantSNPNames = []
	lastPval = 0
	featuresPosList = []
	for l in csvReader:
		snpName, pVal = l[0], float(l[4])
		if (pVal < lastPval): raise Exception('P-values are not sorted in descending order: ' + str(pVal) + ">" + str(lastPval))
		lastPval = pVal
		if (pVal > pthresh): break		
		if (snpName not in snpNameToNumDict): continue							
		significantSNPNames.append(snpName)
		if (mindist == 0):
			significantSNPs.append(snpNameToNumDict[snpName])
			print 'Using SNP', snpName, 'with p<%0.2e'%pVal, 'as a fixed effect'
		else:
			posArr = bed.pos[snpNameToNumDict[snpName]]
			chrom, pos = posArr[0], int(posArr[2])				
			addSNP = True
			for (c,p) in featuresPosList:
				if (chrom == c and abs(pos-p) < mindist):
					addSNP = False
					break
			if addSNP:
				significantSNPs.append(snpNameToNumDict[snpName])
				featuresPosList.append((chrom, pos))
				print 'Using SNP', snpName, '('+str(int(chrom))+':'+str(pos)+') with p<%0.2e'%pVal, 'as a fixed effect'
	f.close()

	snpCovarsMat = bed.val[:, significantSNPs]
	return snpCovarsMat
	
	
	
def getExcludedChromosome(bfile, chrom):
	bed = Bed(bfile)	
	indsToKeep = (bed.pos[:,0] != chrom)
	bed = bed[:, indsToKeep]	
	return bed.read().standardize()
	
def getChromosome(bfile, chrom):
	bed = Bed(bfile)
	indsToKeep = (bed.pos[:,0] == chrom)
	bed = bed[:, indsToKeep]	
	return bed.read().standardize()
	

def _fixupBedAndPheno(bed, pheno, missingPhenotype='-9'):
	bed = _fixupBed(bed)
	bed, pheno = _fixup_pheno(pheno, bed, missingPhenotype)
	return bed, pheno
	
def _fixupBed(bed):
	if isinstance(bed, str):
		return Bed(bed).read().standardize()
	else: return bed

def _fixup_pheno(pheno, bed=None, missingPhenotype='-9'):
	if (isinstance(pheno, str)):
		if (bed is not None):
			bed, pheno = loadPheno(bed, pheno, missingPhenotype, keepDict=True)
			return bed, pheno
		else:
			phenoDict = phenoUtils.loadOnePhen(pheno, missing=missingPhenotype, vectorize=True)
			return phenoDict
	else:
		if (bed is not None): return bed, pheno			
		else: return pheno

def linreg(bed, pheno):

	#Extract snps and phenotype
	bed, pheno = _fixupBedAndPheno(bed, pheno)	
	if isinstance(pheno, dict):	phe = pheno['vals']	
	else: phe = pheno		
	if (len(phe.shape)==2):
		if (phe.shape[1]==1): phe=phe[:,0]
		else: raise Exception('More than one phenotype found')	

	#Normalize y. We assume X is already normalized.
	y = phe - phe.mean(); y /= y.std()

	#Compute p-values
	Xy = bed.val.T.dot(y) / y.shape[0]
	Xy[Xy>1.0] = 1.0
	Xy[Xy<-1.0] = -1.0
	df = y.shape[0]-2
	TINY = 1.0e-20
	t = Xy * np.sqrt(df / ((1.0-Xy+TINY) * (1.0+Xy+TINY)))
	pValT = stats.t.sf(np.abs(t), df)*2	
	
	#Create pandas data frame
	items = [
		('SNP', bed.sid),
		('Chr', bed.pos[:,0]), 
		('GenDist', bed.pos[:,1]),
		('ChrPos', bed.pos[:,2]), 
		('PValue', pValT),                
	]
	frame = pd.DataFrame.from_items(items)	
	frame.sort("PValue", inplace=True)
	frame.index = np.arange(len(frame))	
	return frame
	
def powerPlot(df, causalSNPs, title=''):
	import pylab
	causalSNPs = set(causalSNPs)
	csnpPvals = df[df['SNP'].isin(causalSNPs)]["PValue"]	
	pvalPoints = np.logspace(-6, -2, num=1000)
	power = [np.mean(csnpPvals < p ) for p in list(pvalPoints)]
	pylab.plot(-np.log10(pvalPoints), power)
	pylab.xlabel("-log10(Significance Threshold)")
	pylab.ylabel("Power")
	pylab.title(title)
	
	
def computeCovar(bed, shrinkMethod, fitIndividuals):
	eigen = dict([])

	if (shrinkMethod in ['lw', 'oas', 'l1', 'cv']):
		import sklearn.covariance as cov
		t0 = time.time()
		print 'Estimating shrunk covariance using', shrinkMethod, 'estimator...'
				
		if (shrinkMethod == 'lw'): covEstimator = cov.LedoitWolf(assume_centered=True, block_size = 5*bed.val.shape[0])
		elif (shrinkMethod == 'oas'): covEstimator = cov.OAS(assume_centered=True)
		elif (shrinkMethod == 'l1'): covEstimator = cov.GraphLassoCV(assume_centered=True, verbose=True)
		elif (shrinkMethod == 'cv'):
			shrunkEstimator = cov.ShrunkCovariance(assume_centered=True)
			param_grid = {'shrinkage': [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]}			
			covEstimator = sklearn.grid_search.GridSearchCV(shrunkEstimator, param_grid)		
		else: raise Exception('unknown covariance regularizer')
		
		covEstimator.fit(bed.val[fitIndividuals, :].T)
		if (shrinkMethod == 'l1'):
			alpha = covEstimator.alpha_
			print 'l1 alpha chosen:', alpha
			covEstimator2 = cov.GraphLasso(alpha=alpha, assume_centered=True, verbose=True)
		else:
			if (shrinkMethod == 'cv'): shrinkEstimator = clf.best_params_['shrinkage']
			else: shrinkEstimator = covEstimator.shrinkage_
			print 'shrinkage estimator:', shrinkEstimator
			covEstimator2 = cov.ShrunkCovariance(shrinkage=shrinkEstimator, assume_centered=True)
		covEstimator2.fit(bed.val.T)
		XXT = covEstimator2.covariance_ * bed.val.shape[1]
		print 'Done in %0.2f'%(time.time()-t0), 'seconds'
			
	else:
		print 'Computing kinship matrix...'	
		t0 = time.time()
		XXT = symmetrize(blas.dsyrk(1.0, bed.val, lower=1))
		print 'Done in %0.2f'%(time.time()-t0), 'seconds'		
		try: shrinkParam = float(shrinkMethod)
		except: shrinkParam = -1
		if (shrinkMethod == 'mylw'):
			XXT_fit = XXT[np.ix_(fitIndividuals, fitIndividuals)]
			sE2R = (np.sum(XXT_fit**2) - np.sum(np.diag(XXT_fit)**2)) / (bed.val.shape[1]**2)
			#temp = (bed.val**2).dot((bed.val.T)**2)
			temp = symmetrize(blas.dsyrk(1.0, bed.val[fitIndividuals, :]**2, lower=1))
			sER2 = (temp.sum() - np.diag(temp).sum()) / bed.val.shape[1]
			shrinkParam = (sER2 - sE2R) / (sE2R * (bed.val.shape[1]-1))		
		if (shrinkParam > 0):
			print 'shrinkage estimator:', 1-shrinkParam
			XXT = (1-shrinkParam)*XXT + bed.val.shape[1]*shrinkParam*np.eye(XXT.shape[0])
	
	return XXT


	
	
def standardize(X, method, optionsDict):
	fitIndividuals = np.ones(X.shape[0], dtype=np.bool)
	if (method == 'frq'):
		empMean = X.mean(axis=0) / 2.0
		X[:, empMean>0.5] = 2 - X[:, empMean>0.5]	
		print 'regularizng SNPs according to frq file...'
		frqFile = (optionsDict['bfilesim']+'.frq' if (optionsDict['frq'] is None) else optionsDict['frq'])
		mafs = np.loadtxt(frqFile, usecols=[1,2]).mean(axis=1)
		snpsMean = 2*mafs
		snpsStd = np.sqrt(2*mafs*(1-mafs))	
	elif (method == 'related'):
		if (optionsDict['related'] is None): raise Exception('related file not supplied')
		print 'regularizng SNPs according to non-related individuals...'
		relLines = np.loadtxt(optionsDict['related'], usecols=[2])	
		keepArr = (relLines != 1)
		print 'Excluding', np.sum(~keepArr), 'from the covariance matrix standardization'
		snpsMean = X[keepArr, :].mean(axis=0)
		snpsStd = X[keepArr, :].std(axis=0)
		fitIndividuals = keepArr
	elif (method == 'controls'):
		phe = optionsDict['pheno']
		pheThreshold = phe.mean()
		controls = (phe<pheThreshold)		
		print 'regularizng SNPs according to controls...'
		snpsMean = X[controls, :].mean(axis=0)
		snpsStd = X[controls, :].std(axis=0)
		fitIndividuals = controls
	elif (method is None):
		snpsMean = X.mean(axis=0)
		snpsStd = X.std(axis=0)
	else:
		raise Exception('unknown SNP standardization option: ' + method)

	X -= snpsMean, 
	X /= snpsStd
	return X, fitIndividuals
