P
----------------

LEAP is a program for liability estimation in ascertained case-control studies, written in the Python language.
LEAP estimates liabilities, that can then be treated as phenotypes in a GWAS context, which can greatly increase power.


------------------
Installation
------------------
The easiest way to install LEAP is via pip, by typing the command:
```shell
pip install --user leap_gwas
```

Typically, the LEAP scripts will be installed at:
```
~/.local/lib/python2.7/site-packages/LEAP/
```

LEAP is particularly easy to install using the [Anaconda Python distribution](https://store.continuum.io/cshop/anaconda). The [numerically optimized version](http://continuum.io/blog/mkl-optimizations) of Anaconda can speed LEAP up by several orders of magnitude.

Alternatively (if Anaconda can't be installed), for very fast performance it is recommended to have an optimized version of Numpy/Scipy [installed on your system](http://www.scipy.org/scipylib/building), using optimized numerical libraries such as [OpenBLAS](http://www.openblas.net) or [Intel MKL](https://software.intel.com/en-us/intel-mkl) (see [Compilation instructions for scipy with Intel MKL)](https://software.intel.com/en-us/articles/numpyscipy-with-intel-mkl). 

If you want to install LEAP manually you need the following dependencies:
* Python 2.7
* [Numpy](http://www.numpy.org/) and [Scipy](http://www.scipy.org/)
* [Scikits-learn](http://scikit-learn.org/stable/)
* The Python [FaST-LMM package](https://github.com/MicrosoftGenomics/FaST-LMM).

Please make sure these are installed prior to using LEAP.
 
 
------------------
Usage instructions
----------------------
LEAP is invoked though a series of Python scripts, as detailed below.
The script leapUtils.sh runs the full LEAP pipeline on a small example dataset, and can be used for reference.
 
Generally, LEAP uses the same file formats as [FaST-LMM](https://github.com/MicrosoftGenomics/FaST-LMM).
Namely, input files are in [binary Plink format](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed).
When there is a contradiction between file formats used by Plink and by FastLMM, LEAP uses the convention adopted by FastLMM.
Explanations about the parameters used by all the scripts can be seen by typing
```
python <script_name> --help
```

The command-line options for LEAP largely follow the options of the [C++ version of 
FaST-LMM](http://research.microsoft.com/en-us/projects/fastlmm/).
 
 
### The LEAP pipeline
**1) (optional): Find related individuals to be removed:**
```
python findRelated.py --bfilesim <Plink base file> --out <output file>
```
 This script creates a file marking the individuals that need to be removed to eliminate relatedness.
 
**2) Compute heritability using the method of [Golan et al.](http://www.pnas.org/content/111/49/E5272.long):**
```
python calc_h2.py --bfilesim <Plink base file> --prev <prevalence> --pheno <phenotype file> --h2coeff 1.0 [--extractSim <SNPs used for heritability estimation>  --related <relatedness file> --h2coeff <heritability coefficient>]
```
 This script outputs the heritability estimate. The optional extractSim file is a text file with a list of SNP names (one SNP per line) that will be used for heritability estimation. It is recommended to perform a different heritability and liability estimation for every excluded chromosome, and then testing the SNPs on the excluded chromosome for association with the estimated liabilities. The optional relatedness file should be the output of stage 1, and is used to exclude related individuals from the analysis, which improves analysis results.

**3) Estimate liabilities:**
```
python probit.py --bfilesim <Plink base file> --pheno <phenotype file> --prev <prevalence> --out <output base file> --h2 <heritability> [--extractSim <SNPs used in the liability estimation> --related <relatedness file>]
```
This script creates a file called \<output base file\>.liabs, with estimated liabilities for every individual. The estimated liabilities can be used directly for GWAS by using them as a standard phenotype file. The h2 parameter should be the heritability estimate from stage 2. The extractSim and relatedness file parameters should be the same as in stage 2.

**4) Test for Associations:**
```
python leap_gwas.py --bfilesim <Plink base file for kinship estimation> --bfile <Plink file with tested SNPs> --pheno <estimated liabilities file> --out <output file> --h2 <heritability> [--extractSim <SNPs used in the LMM kinship matrix>  --extract <SNPs to test>]
```
This script performs GWAS with a prespecified heritability level (as computed in stage 2). The pheno parameter is the liabilities file computed in stage 3. The syntax largely follows that of the [C++ version of FaST-LMM](http://research.microsoft.com/en-us/projects/fastlmm/).
Kinship estimation should use the same SNPs used for heritability and liability estimation. The bfile and bfilesim parameters can both point to the same file. In this case, the extract and extractSim parameters should be used to guarantee that kinship estimation doesn't use SNPs on the excluded chromosome, and that all tested SNPs are on the excluded chromosome.
 
 
-----------------
General comments and tips
-------------------------
**1)** Fixed effects can be included in stages 2-4 by adding the flag --covar.
Please type
```
python probit.py --help
```
and
```
python leap_gwas.py --help
```
for instructions. However, we note that under extreme ascertainment, it is recommded to use covariates only in stages 2-3 (see the paper for details).
 
**2)** As described in the main text, it is recommended to perform a different liability estimation for every excluded chromosome, and then testing the SNPs on the excluded chromosome for association with the estimated liabilities. The -extractSim flag is useful for this. Please see the example file leap_pipeline.sh for a usage example.
 
**3)** A complete end-to-end usage example is provided with the LEAP source files, and can be invoked via the script leap_pipeline.sh.
This example estimates liabilities for a small balanced case-control dataset.
The dataset was simulated with 50% heritability and  0.1% prevalence. It included 500 cases, 500 controls, 499 causal SNPs, 100 unusually differentiated SNPs and 10000 SNPs differentiated with FST=0.01. Causal SNPs are called csnp\<i\>, and unusually differentiated SNPs are called dsnp\<i\>. The original liabilities for this file are available in the file dataset1.phe.liab (but this file is not used by LEAP).
 

-----------------
Contact
---------
For questions and comments, please contact Omer Weissbrod at omerw[at]cs.technion.ac.il


