#!/bin/bash

##申请资源
#SBATCH -J PROSPERV2
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 64

#进入虚拟环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py36

#需要执行的命令
rootdir=/work/home/acd2j8na2s/Work/PRSprofile
Rpwd=/work/home/acd2j8na2s/soft/R/4.0.3/bin/Rscript
cd ${rootdir}
plinkpwd2='/work/home/acd2j8na2s/software/plink2'


cd /work/home/acd2j8na2s/software/PROSPER
for trait in SBP
do
		$Rpwd /work/home/acd2j8na2s/software/PROSPER/scripts/PROSPER.revise.R \
		--PATH_package /work/home/acd2j8na2s/software/PROSPER \
		--PATH_out ${rootdir}/7_single_ancestry_PRS/$trait/PROSPER \
		--FILE_sst ${rootdir}/7_single_ancestry_PRS/$trait/PROSPER/EAS.txt,${rootdir}/7_single_ancestry_PRS/$trait/PROSPER/EUR.txt \
		--pop EAS,EUR \
		--lassosum_param ${rootdir}/7_single_ancestry_PRS/$trait.EASoptimal_param.txt,${rootdir}/7_single_ancestry_PRS/$trait.EURoptimal_param.txt \
		--chrom 1-22 \
		--NCORES 10
done
#Error in { : 
#  task 15 failed - "number of columns of matrices must match (see arg 2)"
#Calls: %dopar% -> <Anonymous>
#Execution halted
#没有15号染色体?why
#ERROR: 15号染色体 in ## Step 2.4. Clean PRSs into a matrix (#variant X #grid_search)
##number of columns of matrices must match (see arg 2)
  #snps_scale[[1/2]] NULL
  #dim(prs)
  #1] 20343    50