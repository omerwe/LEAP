###    LEAP

LEAP is a program for liability estimation in ascertained case-control studies, written in the Python language.
It can estimate liabilities, that can then be treated as phenotypes in a GWAS context, which can greatly increase power.


###    Installation
* The easiest way to install LEAP is via pip, by typing the command:
```shell
pip install --user leap_gwas
```

Typically, the LEAP scripts will be installed at:
```
~/.local/lib/python2.7/site-packages/LEAP/
```

* LEAP is particularly easy to install using the [Anaconda Python distribution](https://store.continuum.io/cshop/anaconda). The [numerically optimized version](http://continuum.io/blog/mkl-optimizations) of Anaconda can speed LEAP up by several orders of magnitude.
* Alternatively (if Anaconda can't be installed), for very fast performance it is recommended to have an optimized version of Numpy/Scipy [installed on your system](http://www.scipy.org/scipylib/building), using optimized numerical libraries such as [OpenBLAS](http://www.openblas.net) or [Intel MKL](https://software.intel.com/en-us/intel-mkl) (see [Compilation instructions for scipy with Intel MKL)](https://software.intel.com/en-us/articles/numpyscipy-with-intel-mkl). 

* If you want to install LEAP from source you need the following dependencies:
* Python 2.7
* Numpy and Scipy
* Scikits-learn
* The pysnptools package (https://github.com/MicrosoftGenomics/PySnpTools).
Please make sure these are installed prior to using LEAP.
 
 

###    Usage instructions
LEAP is invoked though a series of Python scripts, as detailed below.
The script leapUtils.sh runs the full LEAP pipeline on a small example dataset, and can be used for reference.
 
Generally, LEAP uses the same file formats as FastLMM.
Namely, input files are in binary Plink format (http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed).
When there is a contradiction between file formats used by Plink and by FastLMM, LEAP uses the convention adopted by FastLMM.
Explanations about the parameters used by all the scripts can be seen by typing
```
python <script_name> --help
```
 
 
The LEAP pipeline includes the following stages:
-------------------------------------------------
1) (optional): Find related individuals to be removed:
```
python findRelated.py --bfilesim <Plink base file> --out <output file>
```
 This script creates a file marking the individuals that need to be removed to eliminate relatedness
 
2) Compute heritability using the method of Golan and Rosset:
```
python calc_h2.py --bfilesim <Plink base file> --extractSim <SNPs used for heritability estimation> --prev <prevalence> --pheno <phenotype file> --related <relatedness file>
```
 This script outputs the heritability estimate
 
3) Estimate liabilities:
```
python probit.py --bfilesim <Plink base file> --pheno <phenotype file> --prev <prevalence> --extractSim <SNPs used in the heritability estimation> --out <output base file> --related <relatedness file> --h2 <heritability>
```
 This script creates a file called \<output base file\>.liabs, with estimated liabilities for every individual. The estimated liabilities can be used directly for GWAS by using them as a standard phenotype file.

4) Compute GWAS:
```
python leap_gwas.py --bfilesim <Plink base file> --pheno <estimated liabilities file> --extractSim <SNPs used in the LMM kinship matrix> --out <output file> --h2 <heritability> --bfile <Plink file with tested SNPs> --extract <SNPs to test>
```
 This script performs GWAS with a prespecified h2. The syntax largely follows that of FaSTLMM C++ version.

 
 
 
General comments and tips
-------------------------
1. Fixed effects can be included in the thresholds, liability estimation and GWAS stages.
Please type
```
python probit.py --help
```
and
```
leap_gwas.py
```
for instructions.
 
2. As described in the main text, it is recommended to perform a different liability estimation for every excluded chromosome, and then testing the SNPs on the excluded chromosome for association with the estimated liabilities. The -extractSim flag is useful for this. Please see the example file leap_pipeline.sh for a usage example.
 
3. A complete end-to-end usage example is provided with the LEAP source files, and can be invoked via the script leap_pipeline.sh.
This example estimates liabilities for a small balanced case-control dataset.
The dataset was simulated with 50% heritability and  0.1% prevalence. It included 500 cases, 500 controls, 499 causal SNPs, 100 unusually differentiated SNPs and 10000 SNPs differentiated with FST=0.01. Causal SNPs are called csnp\<i\>, and unusually differentiated SNPs are called dsnp\<i\>. The original liabilities for this file are available in the file dataset1.phe.liab (but this file is not used by LEAP).
 


## Contact
For questions and comments, please contact Omer Weissbrod at omerw[at]cs.technion.ac.il

