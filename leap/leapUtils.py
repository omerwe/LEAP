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
	covarsDict = phenoUtils.loadOnePhen(covarFile, vectorize=False)
	checkIntersection(bed, covarsDict, 'covariates', checkSuperSet=True)
	_, covarsDict = pstutil.intersect_apply([bed, covarsDict])
	covar = covarsDict['vals']
	covar -= np.mean(covar, axis=0)
	covar /= np.std(covar, axis=0)	
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
	
	
	