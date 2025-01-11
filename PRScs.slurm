#!/bin/bash

##申请资源
#SBATCH -J PRScs
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
  for db in BBJ UKBB
  do
  if [ $db == BBJ ]
  then
    REFDIR=$rootdir/PRScs/ldblk_1kg_eas
    if [ $trait == BMI ]
    then
    N_gwas=158284
    elif [ $trait == SBP ]
    then
    N_gwas=136597
    else
    N_gwas=136615
    fi
  else
    REFDIR=$rootdir/PRScs/ldblk_1kg_eur
    if [ $trait == BMI ]
    then
    N_gwas=457824
    elif [ $trait == SBP ]
    then
    N_gwas=422771
    else
    N_gwas=422771
    fi
  fi
  echo $REFDIR
  echo $N_gwas
  python ${rootdir}/PRScs/PRScs-master/PRScs.py \
  --ref_dir=${REFDIR} \
  --bim_prefix=${rootdir}/5_allGWASandData/$trait/plink/STS2030 \
  --sst_file=$rootdir/5_allGWASandData/$trait/$db/$trait.gwas.QC.CAS.PRScs.txt \
  --n_gwas=$N_gwas \
  --out_dir=$rootdir/7_single_ancestry_PRS/$trait/PRScs/$db
  done
done
