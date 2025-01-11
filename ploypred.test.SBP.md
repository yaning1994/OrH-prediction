# 1. Genome-wide fine-mapping using PolyFun
This step performs fine-mapping in each locus in the genome, using PolyFun. Please read the PolyFun documentation for details. Here is a brief end-to-end example that you can try out:
```sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate polyfun
mkdir /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun
cd /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun
```
```sh
#created munged sumstats
cat /work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/SBP/UKBB/SBP.gwas.QC.CAS | head
#SNP     CHR     BP      EA      OA      EAF     Beta    se      P       N

python /work/home/acd2j8na2s/Work/polyfun/polyfun/munge_polyfun_sumstats.py \
    --sumstats /work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/SBP/UKBB/SBP.gwas.QC.CAS \
    --out /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/bolt.sumstats.munged.parquet
#ValueError: cannot both specify --n and have an N column in the sumstats file

#in python read 
#df = pd.read_parquet('/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/bolt.sumstats.munged.parquet')
#SNP  CHR      BP A1 A2       MAF       N         Z
# N= 457824
```
```sh
#Using precomputed prior causal probabilities based on a meta-analysis of 15 UK Biobank traits
python /work/home/acd2j8na2s/Work/polyfun/polyfun/extract_snpvar.py \
--sumstats /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/bolt.sumstats.munged.parquet \
--out /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/snps_with_var.gz
#sbatch /work/home/acd2j8na2s/Work/PRSprofile/0_script/polyfunSNPVAR.slurm
```

```sh
cd /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun
#create fine-mapping jobs
mkdir /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step1
python /work/home/acd2j8na2s/Work/polyfun/polyfun/create_finemapper_jobs.py \
    --sumstats /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/snps_with_var.gz \
    --n 457824 \
    --method susie \
    --max-num-causal 5 \
    --out  /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step1/polyfun_output \
    --jobs-file /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/jobs.txt \
    --python3 /work/home/acd2j8na2s/miniconda3/envs/polyfun/bin/python
```
polyfunjobs.slurm中replace --ld https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/ as --ld /work/home/acd2j8na2s/Work/polyfun/LD_chr_SE/
否则#urllib.error.URLError: <urlopen error [Errno 110] Connection timed out>,而且每次不完全一致
```sh
cd LD_chr_SE
#total 5408 file
ls /work/home/acd2j8na2s/Work/polyfun/LD_chr_SE | wc -w
cd ..
```
下载总是失败，使用Chrono是谷歌浏览器里的扩展程序批量下载到本地后上传

```sh
#run each and every fine-mapping job (could take several hours or more)
sbatch /work/home/acd2j8na2s/Work/PRSprofile/0_script/polyfunjobs.slurm
#jobs id 4586547
```

```sh
#aggregate all of the results
python /work/home/acd2j8na2s/Work/polyfun/polyfun/aggregate_finemapper_results.py \
    --out-prefix /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step1/polyfun_output \
    --sumstats /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/snps_with_var.gz \
    --out /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step1/polyfun_output.agg.txt.gz \
    --adjust-beta-freq \
    --allow-missing-jobs
#The flag --adjust-beta-freq adjusts the posterior effect size estimates to be on a per-causal-allele scale rather than a per-standardized-genotype scale. 
#sbatch /work/home/acd2j8na2s/Work/PRSprofile/0_script/polyfunjobs.aggregate.slurm
#4674973 ing
```
# 2. Estimating tagging SNP effect sizes using another method
这一步不需要额外做，因为GWAS结果本身就是基于BOLT-LMM处理的
这里只需要把生成bolt.betas.gz文件
```sh
#看一下示例文件的格式
#cat /work/home/acd2j8na2s/Work/polyfun/polyfun/polypred_example/bolt.betas.gz | zcat | head 
#SNP     CHR     BP      GENPOS  ALLELE1 ALLELE0 BETA
mkdir /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step2
```
用下边的文件试试
```py
import pandas as pd 
import numpy as np 
#10      rs9419461       124767  T       C       -6.533581e-07
df = pd.read_csv('/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/PRScs/UKBB/pst_eff_a1_b0.5_phiauto.txt',sep='\t',header=None)
df.columns=['CHR','SNP','BP','EA','OA','Beta']
df['GENPOS']=0
df = df[['SNP','CHR','BP','GENPOS','EA','OA','Beta']]
df.columns=['SNP','CHR','BP','GENPOS','ALLELE1','ALLELE0','BETA']
df.to_csv('/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step2/uk.PRScs.betas',sep='\t',index=False)

df = pd.read_csv('/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/PRScs/BBJ/pst_eff_a1_b0.5_phiauto.txt',sep='\t',header=None)
df.columns=['CHR','SNP','BP','EA','OA','Beta']
df['GENPOS']=0
df = df[['SNP','CHR','BP','GENPOS','EA','OA','Beta']]
df.columns=['SNP','CHR','BP','GENPOS','ALLELE1','ALLELE0','BETA']
df.to_csv('/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step2/bbj.PRScs.betas',sep='\t',index=False)
```
# 3. Linearly combining the effect sizes of PolyFun and the other method in small training set
```sh
cat /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step2/uk.PRScs.betas | head 
#SNP     CHR     BP      GENPOS  ALLELE1 ALLELE0 BETA

mkdir /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step3
```

source ~/miniconda3/etc/profile.d/conda.sh
conda activate py36
#需要执行的命令
```R
#/work/home/acd2j8na2s/soft/R/4.0.3/bin/R
#phenofile: FID    IID PHENO
library(bigsnpr)
library(bigreadr)
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(boot)


rootdir <- '/work/home/acd2j8na2s/Work/PRSprofile'
db <- 'UKBB'
trait <- 'SBP'
restdir <- '/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step3'
genfile <- paste(rootdir,'5_allGWASandData',trait,'plink',sep = "/")
plinkfile <- paste(genfile,'ALL',sep = "/")
pcafile <- paste(genfile,'ALL.eigenvec',sep = "/")
phenofilefile <- paste(rootdir,'7_single_ancestry_PRS/ALL.phe.txt',sep = "/")
if(trait=='SBP') {
    affection <- 'sbp'
}
if(trait=='SBP') {
    affection <- 'SBP'
}
if(trait=='SBP') {
    affection <- 'sbp'
}

#read in phenotype and covariates
phenotype <- fread(phenofilefile,sep='\t')
phenotype$FID <- phenotype$id
phenotype$IID <- phenotype$id
phenotype <- phenotype %>% select("FID", "IID", 'age','sex',affection) 
#data %>% select(all_of(affection))
names(phenotype) <- c("FID", "IID", 'age','sex',"target")

pcs <- fread(pcafile)
# rename columns
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
# generate required table
pheno <- merge(phenotype, pcs)
#Y~age+sex+6pc
formulaagj <- paste("target","~age+sex+PC1+PC2+PC3+PC4+PC5+PC6")
fit <- lm(formulaagj, dat=pheno)
resid <- residuals(fit)
pheno$target <- resid
#read in famfile
fam.order <- fread(paste(plinkfile,".fam",sep=''))
fam.order <- fam.order %>% select("V1", "V2") 
# We assume the fam order is the same across different chromosomes
# Rename fam order
setnames(fam.order,
        c("V1", "V2"),
        c("FID", "IID"))
y <- pheno[fam.order, on = c("FID", "IID")]
sts2030fam <- read.table(
  file = "/work/home/acd2j8na2s/Work/PRSprofile/1_Data/STS2030.fam", 
  sep = " ", header = F)
sts1071fam <- read.table(
  file = "/work/home/acd2j8na2s/Work/PRSprofile/1_Data/STS1071.fam", 
  sep = " ", header = F)
caspmifam <- read.table(
  file = "/work/home/acd2j8na2s/Work/PRSprofile/1_Data/CASPMI.fam", 
  sep = " ", header = F)

ind.val  <- which(fam.order$IID %in% c(sts2030fam$V1))
ind.test  <- which(fam.order$IID %in% c(caspmifam$V1))
ind.extend  <- which(fam.order$IID %in% c(sts1071fam$V1))
#Phenotypes for tuning the target PRS 
y <- rename(y, c("PHENO"="target"))
y_tune <- y[ind.val]%>% select("FID", "IID","PHENO") 
#Phenotypes for validation the target PRS 
y_vad <- y[ind.test]%>% select("FID", "IID","PHENO") 
#Phenotypes for test the target PRS 
y_test <- y[ind.extend]%>% select("FID", "IID","PHENO") 
##phenofile: FID    IID PHENO
write.table (y_tune, file =paste(restdir,'y_tune.txt',sep='/'), sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
write.table (y_vad, file =paste(restdir,'y_vad.txt',sep='/'), sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
write.table (y_test, file =paste(restdir,'y_test.txt',sep='/'), sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
```
```sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate polyfun
cd /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun
#genofile: 1000G.subset.$chr.bed/bim/fam
plinkpwd=/work/home/acd2j8na2s/software/plink
mkdir /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step3/plink 
for i in {1..22}
do
$plinkpwd --bfile /work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/SBP/plink/STS2030 \
--chr $i \
 --make-bed -out \
 /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step3/plink/STS2030.$i
done
python /work/home/acd2j8na2s/Work/polyfun/polyfun/polypred.py \
    --combine-betas \
    --betas /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step2/uk.PRScs.betas,/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step2/bbj.PRScs.betas,/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step1/polyfun_output.agg.txt.gz \
    --pheno /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step3/y_tune.txt \
    --output-prefix /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step3/polypredP \
    --plink-exe /work/home/acd2j8na2s/software/plink \
    /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step3/plink/STS2030.*.bed 

#sbatch /work/home/acd2j8na2s/Work/PRSprofile/0_script/polyfunjobs.minbeta.slurm
#4525320
#[INFO]  Writing mixing weights to /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step3/polypred.mixweights
#[INFO]  Loading betas file /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step2/uk.bolt.betas...
#[INFO]  done in 4.18 seconds
#[INFO]  Loading betas file /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step2/bbj.bolt.betas...
#[INFO]  done in 4.26 seconds
#[INFO]  Loading betas file /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step1/polyfun_output.agg.txt.gz...
#[INFO]  done in 12.72 seconds
#[INFO]  Saving weighted betas to /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/SBP/polyfun/step3/polypred.betas
```

# 4. Computing PRS in the target population using the combined SNP effect sizes
```sh
#进入虚拟环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py36

#需要执行的命令
rootdir=/work/home/acd2j8na2s/Work/PRSprofile
cd ${rootdir}
plinkpwd=/work/home/acd2j8na2s/software/plink
Rpwd=/work/home/acd2j8na2s/soft/R/4.0.3/bin/Rscript

for trait in SBP
do
    $Rpwd $rootdir/0_script/PRSinRforPolypred.R $trait
done
```

