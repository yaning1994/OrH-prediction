#!/bin/bash

##申请资源
#SBATCH -J PRScsx
#SBATCH -p xahctest
#SBATCH -N 1
#SBATCH -n 20

#进入虚拟环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py36

#需要执行的命令
rootdir=/work/home/acd2j8na2s/Work/PRSprofile
cd ${rootdir}
plinkpwd=/work/home/acd2j8na2s/software/plink
Rpwd=/work/home/acd2j8na2s/soft/R/4.0.3/bin/Rscript
for trait in BMI DBP SBP
do
    REFDIR=$rootdir/PRScs
    if [ $trait == BMI ]
    then
    N_gwas1=158284
    N_gwas2=457824
    elif [ $trait == SBP ]
    then
    N_gwas1=136597
    N_gwas2=422771
    else
    N_gwas1=136615
    N_gwas2=422771
    fi
  echo $REFDIR
  python $rootdir/0_script/PRScsx-master/PRScsx.py \
  --ref_dir=${REFDIR} \
  --bim_prefix=${rootdir}/5_allGWASandData/$trait/plink/STS2030 \
  --sst_file=${rootdir}/5_allGWASandData/$trait/BBJ/$trait.gwas.QC.CAS.PRScs.txt,${rootdir}/5_allGWASandData/$trait/UKBB/$trait.gwas.QC.CAS.PRScs.txt \
  --n_gwas=${N_gwas1},${N_gwas2} \
  --pop=EAS,EUR \
  --meta=True \
  --out_dir=${rootdir}/7_single_ancestry_PRS/$trait/PRScsx \
  --out_name=PRScsx
done