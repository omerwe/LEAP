#!/bin/bash
mkdir -p results


############ Liability Estimation Part ###############

#Find related individuals to remove
python ../findRelated.py --bfilesim dataset1/dataset1 --out results/dataset1.related

#Compute eigendecompositions for each left-out chromosome
seq 1 10 | xargs -i --max-procs=2 bash -c "python ../eigenDecompose.py --bfilesim dataset1/dataset1 --extractSim dataset1/extracts/nochr{}_extract.txt --pheno dataset1/dataset1.phe --out results/dataset1_nochr{}.npz"

#Compute heritability estimates for every excluded chromosome
seq 1 10 | xargs -i --max-procs=2 bash -c "python ../calc_h2.py --bfilesim dataset1/dataset1 --extractSim dataset1/extracts/nochr{}_extract.txt --prev 0.001 --numRemovePCs 10 --pheno dataset1/dataset1.phe --related results/dataset1.related --h2coeff 1.0 --eigen results/dataset1_nochr{}.npz | tail -1 | cut -d\" \" -f2 > results/dataset1_nochr{}.h2"

# #Estimate liabilities
seq 1 10 | xargs -i --max-procs=2 bash -c "python ../probit.py --bfilesim dataset1/dataset1 --pheno dataset1/dataset1.phe --prev 0.001 --extractSim dataset1/extracts/nochr{}_extract.txt --out results/dataset1_nochr{} --related results/dataset1.related --h2 \`cat results/dataset1_nochr{}.h2\` --hess 0  --eigen results/dataset1_nochr{}.npz"

# GWAS for each pseudo-chromosome
seq 1 10 | xargs -i --max-procs=2 bash -c "python ../leap_gwas.py --bfilesim dataset1/dataset1 --bfile dataset1/dataset1 --pheno results/dataset1_nochr{}.liabs --extractSim dataset1/extracts/nochr{}_extract.txt --out results/dataset1_nochr{}.gwas.out.txt --h2 \`cat results/dataset1_nochr{}.h2\` --extract dataset1/extracts/chr{}_extract.txt --eigen results/dataset1_nochr{}.npz"

#Merge the results
head -1 results/dataset1_nochr1.gwas.out.txt > results/dataset1.gwas.out.txt
sort -gk5 results/dataset1_nochr*.gwas.out.txt | grep -v PValue >> results/dataset1.gwas.out.txt


