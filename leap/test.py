import unittest
import logging
import os
import sys
import pandas as pd
import numpy as np


class TestLeap(unittest.TestCase):

	tempout_dir = "tempout"
	gold_dir = "../results_gold"

	@classmethod
	def setUpClass(self):
		from fastlmm.util.util import create_directory_if_necessary
		create_directory_if_necessary(self.tempout_dir, isfile=False)
		self.pythonpath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__))))
		self.bedbase = os.path.join(self.pythonpath, '../dataset1/dataset1')		
		self.phen_fn = os.path.join(self.pythonpath, '../dataset1/dataset1.phe')
		
		#Create eigendecompositions
		logging.info("Creating eigendecomposition files")
		import eigenDecompose
		for i in xrange(1,11):
			output_file = os.path.abspath(os.path.join(self.tempout_dir, 'dataset1_nochr{}.npz'.format(i)))
			extractSim = '../dataset1/extracts/nochr{0}_extract.txt'.format(i)
			eigenDecompose.eigenDecomposeMain(self.bedbase, output_file, extractSim=extractSim)
		
		
		
	# # def test_eigenDecompose(self):
		# # logging.info("Leap test_eigenDecompose")
		# # import eigenDecompose
		
		# # for i in xrange(1,11):
			# # output_file = os.path.abspath(os.path.join(self.tempout_dir, 'dataset1_nochr{}.npz'.format(i)))
			# # ref_file = os.path.abspath(os.path.join(self.gold_dir, 'dataset1_nochr{}.npz'.format(i)))
			# # extractSim = '../dataset1/extracts/nochr{0}_extract.txt'.format(i)
			# # eigenDecompose.eigenDecomposeMain(self.bedbase, output_file, extractSim=extractSim)
			
			# # npz1 = np.load(ref_file)
			# # npz2 = np.load(output_file)
			# # assert np.max(np.abs((npz1['arr_1']-npz2['arr_1']))) < 1e-8, 'Non-matching eigenvalues found'
			
			# # U1 = npz1['arr_0']
			# # U2 = npz2['arr_0']
			# # U1 *= [(1 if U1[0,i]>0 else -1) for i in xrange(U1.shape[1])]	
			# # U2 *= [(1 if U2[0,i]>0 else -1) for i in xrange(U2.shape[1])]
			# # assert np.max(np.abs((U1-U2))) < 1e-8, 'Non-matching eigenvectors found'
		
		
	def test_findrelated(self):
		logging.info("Leap test_findrelated")
		import findRelated
		ref_file = os.path.abspath(os.path.join(self.gold_dir, 'dataset1.related'))
		output_file = os.path.abspath(os.path.join(self.tempout_dir, 'dataset1.related'))
		findRelated.findRelatedMain(self.bedbase, output_file, self.phen_fn)
		self.compare_pheno(ref_file, output_file)
		
	def test_calch2(self):
		logging.info("Leap test_calch2")
		import calc_h2		
		related_file = os.path.abspath(os.path.join(self.gold_dir, 'dataset1.related'))
		
		for i in xrange(1,11):
			h2_file = os.path.abspath(os.path.join(self.gold_dir, 'dataset1_nochr{0}.h2'.format(i)))
			eigen_file = os.path.abspath(os.path.join(self.tempout_dir, 'dataset1_nochr{}.npz'.format(i)))
			extractSim = '../dataset1/extracts/nochr{0}_extract.txt'.format(i)
			expected_h2 = np.loadtxt(h2_file, usecols=[0])
			h2 = calc_h2.calcH2Main(self.bedbase, self.phen_fn, 0.001, relatedFile=related_file, extractSim=extractSim, h2coeff=1.0, eigen=eigen_file)
			assert np.abs(h2-expected_h2)<1e-5, 'Incorrect heritability estimated'
			
			
	def test_probit(self):
		logging.info("Leap test_probit")
		import probit		
		related_file = os.path.abspath(os.path.join(self.gold_dir, 'dataset1.related'))
		
		for i in xrange(1,11):
			h2_file = os.path.abspath(os.path.join(self.gold_dir, 'dataset1_nochr{0}.h2'.format(i)))
			h2 = np.loadtxt(h2_file, usecols=[0])
			eigen_file = os.path.abspath(os.path.join(self.tempout_dir, 'dataset1_nochr{}.npz'.format(i)))
			ref_file = os.path.abspath(os.path.join(self.gold_dir, 'dataset1_nochr{}.liabs'.format(i)))
			extractSim = '../dataset1/extracts/nochr{0}_extract.txt'.format(i)
			output_file = os.path.abspath(os.path.join(self.tempout_dir, 'dataset1_nochr{}'.format(i)))
			probit.probitMain(self.bedbase, self.phen_fn, h2=h2, prev=0.001, extractSim=extractSim, outFile=output_file, relatedFile=related_file, hess=0, eigen=eigen_file)
			self.compare_pheno(ref_file, output_file+'.liabs')
			
						
	def test_gwas(self):
		logging.info("Leap test_gwas")
		import leap_gwas
		
		for i in xrange(1,11):
			h2_file = os.path.abspath(os.path.join(self.gold_dir, 'dataset1_nochr{0}.h2'.format(i)))
			h2 = np.loadtxt(h2_file, usecols=[0])
			ref_file = os.path.abspath(os.path.join(self.gold_dir, 'dataset1_nochr{}.gwas.out.txt'.format(i)))
			extractSim = '../dataset1/extracts/nochr{0}_extract.txt'.format(i)
			extract = '../dataset1/extracts/chr{0}_extract.txt'.format(i)
			output_file = os.path.abspath(os.path.join(self.tempout_dir, 'dataset1_nochr{}.gwas.out.txt'.format(i)))
			liab_file = os.path.abspath(os.path.join(self.gold_dir, 'dataset1_nochr{}.liabs'.format(i)))
			eigen_file = os.path.abspath(os.path.join(self.tempout_dir, 'dataset1_nochr{}.npz'.format(i)))
			leap_gwas.gwasMain(self.bedbase, self.bedbase, liab_file, h2, output_file, extractSim=extractSim, extract=extract, eigen=eigen_file)
			self.compare_gwas(ref_file, output_file)
		
		
	def compare_pheno(self, file1, file2, delimiter=' ', header=None):
		frame1 = pd.read_csv(file1, delimiter=delimiter, header=header)
		frame2 = pd.read_csv(file2, delimiter=delimiter, header=header)
		frame1.columns = ['fid', 'iid', 'pheno']
		frame2.columns = ['fid', 'iid', 'pheno']		
		assert len(frame1) == len(frame2), '# of lines differs from file "{0}"'.format(file2)
		for row_i, row1 in frame1.iterrows():			
			row2 = frame2[(frame2['fid']==row1['fid']) & (frame2['iid']==row1['iid'])].iloc[0]
			for k in row1.keys():
				assert row1[k]==row2[k], "different data for individual '{0}' '{1}'".format(row1['fid'], row1['iid'])
				
				
	def compare_gwas(self, file1, file2):
	
		ref_frame = pd.read_csv(file1, delimiter='\s', header=0, comment=None)
		frame = pd.read_csv(file2, delimiter='\s', header=0, comment=None)
		assert len(ref_frame) == len(frame), "# of pairs differs from file '{0}'".format(file1)
		for _, row in ref_frame.iterrows():			
			sid = row['SNP']
			pvalue = (frame[frame['SNP'] == sid].iloc[0])['PValue']
			assert abs(row['PValue'] - pvalue) < 1e-5, "snp {0} differs too much from file '{1}'".format(sid,file1)


        				
				
				
	

def getTestSuite():    

	suite1 = unittest.TestLoader().loadTestsFromTestCase(TestLeap)
	return unittest.TestSuite([suite1])

if __name__ == '__main__':
    
	suites = unittest.TestSuite([getTestSuite()])    
	r = unittest.TextTestRunner(failfast=True)
	r.run(suites)

logging.info("done with testing")
