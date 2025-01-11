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
cd /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun
python /work/home/acd2j8na2s/Work/polyfun/polyfun/aggregate_finemapper_results.py \
    --out-prefix /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step1/polyfun_output \
    --sumstats /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/snps_with_var.gz \
    --out /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step1/polyfun_output.agg.txt.gz \
    --adjust-beta-freq \
    --allow-missing-jobs
    