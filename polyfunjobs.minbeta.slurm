#!/bin/bash

##申请资源
#SBATCH -J ployfunjobs
#SBATCH -p 	xahcnormal
#SBATCH -N 1
#SBATCH -n 64

#进入虚拟环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate polyfun


#需要执行的命令
cd /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/BMI/polyfun
python /work/home/acd2j8na2s/Work/polyfun/polyfun/polypred.py \
    --combine-betas \
    --betas /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/BMI/polyfun/step2/uk.bolt.betas,/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/BMI/polyfun/step2/bbj.bolt.betas,/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/BMI/polyfun/step1/polyfun_output.agg.txt.gz \
    --pheno /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/BMI/polyfun/step3/y_tune.txt \
    --output-prefix /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/BMI/polyfun/step3/polypred \
    --plink-exe /work/home/acd2j8na2s/software/plink \
    /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/BMI/polyfun/step3/plink/STS2030.*.bed 
    