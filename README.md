##########################
#    Dependencies:       #
##########################
LEAP requires Python 2.7, and is dependent on the FastLMM Python package.
Please make sure these are installed prior to using LEAP

 
 
################################
#    Usage instructions:       #
################################
LEAP is invoked though a series of Python scripts, as detailed below.
The script leapUtils.sh runs the full LEAP pipeline on a small example dataset, and can be used for reference.
 
Generally, LEAP uses the same file formats as FastLMM.
Namely, input files are in binary Plink format (http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed).
When there is a contradiction between file formats used by Plink and by FastLMM, LEAP uses the convention adopted by FastLMM.
Explanations about the parameters used by all the scripts can be seen by typing "python <script_name> --help".
 
 
The LEAP pipeline includes the following stages:
-------------------------------------------------
1. (optional): Find related individuals to be removed:
python findRelated.py --bfilesim <Plink base file> --out <output file>
 
This script creates a file marking the individuals that need to be removed to eliminate relatedness
 
2. Compute heritability using the method of Golan and Rosset:
python calc_h2.py --bfilesim <Plink base file> --extractSim <SNPs used for heritability estimation> --prev <prevalence> --pheno <phenotype file> --related <relatedness file>
 
This script outputs the heritability estimate
 
3. Estimate liabilities:
python probit.py --bfilesim <Plink base file> --pheno <phenotype file> --prev <prevalence> --extractSim <SNPs used in the heritability estimation> --out <output base file> --related <relatedness file> --h2 <heritability>
 
This script creates a file called <output base file>.liabs, with estimated liabilities for every individual. The estimated liabilities can be used directly for GWAS by using them as a standard phenotype file.
These liability files should be used as phenotype files in an LMM-based GWAS, with the parameter logdelta set as log(1/h2 - 1), where h2 is the heritability estimate from stage 2.
 
 
 
General comments and tips
-------------------------
1. When performing GWAS with estimated liabilities, it is recommended to use the value of logdelta corresponding to the estimated liabilities.
For FastLMM C++ version, this is done via the -logdelta flag.
For a given h2 estimate, the corresponding logdelta estimate is: log(1/h2 - 1)
 
2. Fixed effects can be included in the thresholds and liability estimation stages.
Please type "python probit.py --help" for instructions
 
3. As described in the main text, it is recommended to perform a different liability estimation for every excluded chromosome, and then testing the SNPs on the excluded chromosome for association with the estimated liabilities. The "-extractSim" flag is useful for this. Please see the example file "leap_pipeline.sh" for a usage example.
 
4. A complete end-to-end usage example is provided with the LEAP source files, and can be invoked via the script "leap_pipeline.sh".
This example estimates liabilities for a small balanced case-control dataset.
The dataset was simulated with 50% heritability and  0.1% prevalence. It included 500 cases, 500 controls, 499 causal SNPs, 100 unusually differentiated SNPs and 10000 SNPs differentiated with FST=0.01. Causal SNPs are called "csnp<i>", and unusually differentiated SNPs are called "dsnp<i>". The original liabilities for this file are available in the file dataset1.phe.liab (but this file is not used by LEAP).
 
