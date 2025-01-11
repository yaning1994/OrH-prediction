#!/bin/bash

##申请资源
#SBATCH -J jointPRS
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 64

#进入虚拟环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate JointPRS

#需要执行的命令
rootdir=/work/home/acd2j8na2s/Work/PRSprofile
cd ${rootdir}

for trait in BMI DBP SBP
do
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
    mkdir ${rootdir}/7_single_ancestry_PRS/$trait/JointPRS
    JointPRS_path=/work/home/acd2j8na2s/software/JointPRS
    reference_path=${rootdir}/PRScs; type=1KG
    bim_path=${rootdir}/5_allGWASandData/$trait/plink; bim_prefix=STS2030
    outcome_path=${rootdir}/7_single_ancestry_PRS/$trait/JointPRS
    pop1=EAS; pop2=EUR
    r1=1; r2=1
    sst1=${rootdir}/5_allGWASandData/$trait/BBJ/$trait.gwas.QC.CAS.PRScs.txt; sst2=${rootdir}/5_allGWASandData/$trait/UKBB/$trait.gwas.QC.CAS.PRScs.txt
    #JointPRS-auto
    python ${JointPRS_path}/JointPRS.py \
        --ref_dir=${reference_path}/${type} \
        --bim_prefix=${bim_path}/${bim_prefix} \
        --pop=${pop1},${pop2} \
        --rho_cons=${r1},${r2} \
        --sst_file=${sst1},${sst2} \
        --n_gwas=${N_gwas1},${N_gwas2} \
        --out_dir=${outcome_path} \
        --out_name=JointPRS_auto_${pop1}_${pop2}_${r1}${r2}_${type}
    #Code for {1e-06,1e-04,1e-02,1e+00}
    for param_phi in {1e-06,1e-04,1e-02,1e+00}
        do
            python ${JointPRS_path}/JointPRS.py \
            --ref_dir=${reference_path}/${type} \
            --bim_prefix=${bim_path}/${bim_prefix} \
            --pop=${pop1},${pop2} \
            --rho_cons=${r1},${r2} \
            --sst_file=${sst1},${sst2} \
            --n_gwas=${N_gwas1},${N_gwas2} \
            --phi=${param_phi} \
            --out_dir=${outcome_path} \
            --out_name=JointPRS_tune_${pop1}_${pop2}_${r1}${r2}_${type}_${param_phi}
        done
done
