


# genofile:
```sh
cd /liufanGroup/zhangyn/1_Data/2.STS/3.after_imputation_plink
plink --bfile /liufanGroup/zhangyn/1_Data/2.STS/3.after_imputation_plink/imputation.3101.updatemap --remove /liufanGroup/zhangyn/1_Data/2.STS/3.after_imputation_plink/imputation.1071.updatemap.id --make-bed -out /liufanGroup/zhangyn/1_Data/2.STS/3.after_imputation_plink/imputation.2030.updatemap
rootdir=/liufanGroup/zhangyn/2_PRSV3
```
/liufanGroup/public_data/CASPMI/QC/CASPMI
/liufanGroup/zhangyn/1_Data/2.STS/3.after_imputation_plink/imputation.1071.updatemap
/liufanGroup/zhangyn/1_Data/2.STS/3.after_imputation_plink/imputation.2030.updatemap

##########################################################################################################################
# step1.rsid.old & rsid.new
```py
import pandas as pd 
import numpy as np 
from functools import reduce
import math
import sys
#input:
caspmifile = '/liufanGroup/public_data/CASPMI/QC/CASPMI'
sts1071file = '/liufanGroup/zhangyn/1_Data/2.STS/3.after_imputation_plink/imputation.1071.updatemap'
sts2030file = '/liufanGroup/zhangyn/1_Data/2.STS/3.after_imputation_plink/imputation.2030.updatemap'
easfile = '/liufanGroup/public_data/1000G/1000G'
#output:
caspmi2EASidtxt = '/liufanGroup/zhangyn/2_PRSV3/1_Data/CASPMI.991'
sts10712EASidtxt = '/liufanGroup/zhangyn/2_PRSV3/1_Data/imputation.1071'
sts20302EASidtxt = '/liufanGroup/zhangyn/2_PRSV3/1_Data/imputation.2030'
def readbimfile(BIMfile):
	df = pd.read_csv(BIMfile+'.bim',sep='\s+',header=None)
	setcolumns = ['CHR','SNP','cm','BP','A1','A2']
	df.columns = setcolumns
	return(df)


def getTarget2EASid(targetfile,easfile,target2EASidtxt):
	targetbimdf = readbimfile(targetfile)
	easbimdf = readbimfile(easfile)
	# merged on key: chr, pos, a1, a2
	mergedbimdf1 = pd.merge(targetbimdf,easbimdf,on=['CHR','BP','A1','A2'],how='inner')
	mergedbimdf1['SNP'] = mergedbimdf1['SNP_y']
	# merged on key: chr, pos, a2, a1
	mergedbimdf2 = pd.merge(targetbimdf,easbimdf,left_on=['CHR','BP','A1','A2'],right_on=['CHR','BP','A2','A1'],how='inner')
	mergedbimdf2['SNP'] = mergedbimdf2['SNP_y']
	#update variant IDs using command :--update-name <filename> [new ID col. number] [old ID col.] 
	mergedbimdf1 = mergedbimdf1[['SNP_x', 'SNP']]
	mergedbimdf2 = mergedbimdf2[['SNP_x', 'SNP']]
	mergedbimdf = mergedbimdf1.append(mergedbimdf2)
	mergedbimdf = mergedbimdf[['SNP','SNP_x']]
	mergedbimdf.columns=['newid','oldid']
	mergedbimdf = mergedbimdf[['oldid','newid']]
	mergedbimdf = mergedbimdf.drop_duplicates()
	#for sts 'rs34351794' & 'exm2226441' have the same chr,pos,a1,a2
	#so in mergedbimdf, 'rs34351794' duplicate in col'newid', need to drop_duplicates
	mergedbimdf = mergedbimdf.drop_duplicates(subset=['newid'],keep='first')
	mergedbimdf.to_csv(target2EASidtxt,sep='\t',header=None,index=False)


getTarget2EASid(caspmifile,easfile,caspmi2EASidtxt)
getTarget2EASid(sts1071file,easfile,sts10712EASidtxt)
getTarget2EASid(sts2030file,easfile,sts20302EASidtxt)
```
# step2.change rsid as 1kg : plink --update-map
```sh
caspmifile='/liufanGroup/public_data/CASPMI/QC/CASPMI'
sts1071file='/liufanGroup/zhangyn/1_Data/2.STS/3.after_imputation_plink/imputation.1071.updatemap'
sts2030file='/liufanGroup/zhangyn/1_Data/2.STS/3.after_imputation_plink/imputation.2030.updatemap'
caspmi2EASidtxt='/liufanGroup/zhangyn/2_PRSV3/1_Data/CASPMI.991'
sts10712EASidtxt='/liufanGroup/zhangyn/2_PRSV3/1_Data/imputation.1071'
sts20302EASidtxt='/liufanGroup/zhangyn/2_PRSV3/1_Data/imputation.2030'
EASgenofile='/liufanGroup/public_data/1000G/1000G'
plink --bfile ${caspmifile} --update-map ${caspmi2EASidtxt} --update-name --make-bed -out /liufanGroup/zhangyn/2_PRSV3/1_Data/CASPMI.991.updatemap
#Duplicate SNPID such as'rs2722786'. in STS genofile
awk '{print $2}' ${sts1071file}.bim > ${sts1071file}.snp 
awk 'a[$0]++' ${sts1071file}.snp > ${sts1071file}.snp.dul
#1097 SNP duplication, exclude
plink --bfile ${sts1071file} -exclude ${sts1071file}.snp.dul --make-bed -out ${sts1071file}.dul
plink --bfile ${sts2030file} -exclude ${sts1071file}.snp.dul --make-bed -out ${sts2030file}.dul
plink --bfile ${sts1071file}.dul --update-map ${sts10712EASidtxt} --update-name --make-bed -out /liufanGroup/zhangyn/2_PRSV3/1_Data/STS.1071.updatemap
plink --bfile ${sts2030file}.dul --update-map ${sts20302EASidtxt} --update-name --make-bed -out /liufanGroup/zhangyn/2_PRSV3/1_Data/STS.2030.updatemap
```
# step3.只保留CASPMI，STS1071 和STS2030中overloap的SNP(三者SNP数目一样)
```sh

CASPMIQCplink2EASid=$rootdir/1_Data/CASPMI.991.updatemap
STS1071QCplink2EASid=$rootdir/1_Data/STS.1071.updatemap
STS2030QCplink2EASid=$rootdir/1_Data/STS.2030.updatemap
awk '{print $2}' ${CASPMIQCplink2EASid}.bim | uniq > ${rootdir}/1_Data/temp.a
awk '{print $2}' ${STS1071QCplink2EASid}.bim | uniq > ${rootdir}/1_Data/temp.b
awk '{print $2}' ${STS2030QCplink2EASid}.bim | uniq > ${rootdir}/1_Data/temp.c
sort ${rootdir}/1_Data/temp.a ${rootdir}/1_Data/temp.b | uniq -d > ${rootdir}/1_Data/temp.1
sort ${rootdir}/1_Data/temp.1 ${rootdir}/1_Data/temp.c | uniq -d > ${rootdir}/1_Data/cas.sts.eas.snp
wc -l ${rootdir}/1_Data/cas.sts.eas.snp
# 4292543
rm ${rootdir}/1_Data/temp.*
# ---output file
CASPMIQCplink2CasStsEassnp=$rootdir/1_Data/CASPMI
STS1071QCplink2CasStsEassnp=$rootdir/1_Data/STS1071
STS2030QCplink2CasStsEassnp=$rootdir/1_Data/STS2030
# ---plink for extract overloap snp
plink --bfile ${CASPMIQCplink2EASid} \
--extract ${rootdir}/1_Data/cas.sts.eas.snp \
--make-bed \
-out ${CASPMIQCplink2CasStsEassnp}
#4292544 variants and 991 people pass filters and QC.
plink --bfile ${STS1071QCplink2EASid} \
--extract ${rootdir}/1_Data/cas.sts.eas.snp \
--make-bed \
-out ${STS1071QCplink2CasStsEassnp}
#4292543 variants and 1071 people pass filters and QC
plink --bfile ${STS2030QCplink2EASid} \
--extract ${rootdir}/1_Data/cas.sts.eas.snp \
--make-bed \
-out ${STS2030QCplink2CasStsEassnp}
#4292543 variants and 2030 people pass filters and QC.
```
# step4.PCA
```sh
mkdir $rootdir/PCA
plinkfile1=$rootdir/1_Data/CASPMI
plinkfile2=$rootdir/1_Data/STS1071
plinkfile3=$rootdir/1_Data/STS2030
plinkfile4=/liufanGroup/public_data/1000G/1000G
awk '{print $2}' ${plinkfile1}.bim > $rootdir/PCA/snp
plink --bfile $plinkfile4 --extract $rootdir/PCA/snp --make-bed -out $rootdir/PCA/1000G
plink --bfile $rootdir/PCA/1000G --bmerge ${plinkfile1}.bed ${plinkfile1}.bim ${plinkfile1}.fam --make-bed -out $rootdir/PCA/temp1
#Error: 26 variants with 3+ alleles present
plink --bfile $rootdir/PCA/1000G --exclude $rootdir/PCA/temp1-merge.missnp --make-bed -out $rootdir/PCA/temp1.2
plink --bfile $plinkfile1 --exclude $rootdir/PCA/temp1-merge.missnp --make-bed -out $rootdir/PCA/temp1.3
plink --bfile $plinkfile2 --exclude $rootdir/PCA/temp1-merge.missnp --make-bed -out $rootdir/PCA/temp1.4
plink --bfile $plinkfile3 --exclude $rootdir/PCA/temp1-merge.missnp --make-bed -out $rootdir/PCA/temp1.5
plink --bfile $rootdir/PCA/temp1.2 --bmerge $rootdir/PCA/temp1.3.bed $rootdir/PCA/temp1.3.bim $rootdir/PCA/temp1.3.fam --make-bed -out $rootdir/PCA/temp1
plink --bfile $rootdir/PCA/temp1 --bmerge $rootdir/PCA/temp1.4.bed $rootdir/PCA/temp1.4.bim $rootdir/PCA/temp1.4.fam --make-bed -out $rootdir/PCA/PCA
plink --bfile $rootdir/PCA/temp1.4 --exclude $rootdir/PCA/PCA-merge.missnp --make-bed -out $rootdir/PCA/temp1.6
plink --bfile $rootdir/PCA/temp1.5 --exclude $rootdir/PCA/temp1-merge.missnp --make-bed -out $rootdir/PCA/temp1.7
plink --bfile $rootdir/PCA/temp1 --bmerge $rootdir/PCA/temp1.6.bed $rootdir/PCA/temp1.6.bim $rootdir/PCA/temp1.6.fam --make-bed -out $rootdir/PCA/temp1.8
plink --bfile $rootdir/PCA/temp1.8 --bmerge $rootdir/PCA/temp1.7.bed $rootdir/PCA/temp1.7.bim $rootdir/PCA/temp1.7.fam --make-bed -out $rootdir/PCA/PCA
plink --bfile $rootdir/PCA/temp1.8 --exclude $rootdir/PCA/PCA-merge.missnp --make-bed -out $rootdir/PCA/temp1.9
plink --bfile $rootdir/PCA/temp1.7 --exclude $rootdir/PCA/PCA-merge.missnp --make-bed -out $rootdir/PCA/temp1.10
plink --bfile $rootdir/PCA/temp1.10 --bmerge $rootdir/PCA/temp1.9.bed $rootdir/PCA/temp1.9.bim $rootdir/PCA/temp1.9.fam --make-bed -out $rootdir/PCA/PCA
#4292487 variants and 6596 people pass filters and QC
plink --bfile $rootdir/PCA/PCA --pca 4 -out $rootdir/PCA/PCA
python $rootdir/0_script/PCA.py
```
##########################################################################################################################

# step1.phenotype for T1,ST1
```sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py36
rootdir=/work/home/acd2j8na2s/Work/PRSprofile
cd ${rootdir}
cd ${rootdir}/2_Phenotype
python $rootdir/0_script/Phe.V2.py 'STS' 'STS1071'
python $rootdir/0_script/Phe.V2.py 'STS' 'STS2030'
python $rootdir/0_script/Phe.V2.py 'CASPMI' 'CASPMI'
python $rootdir/0_script/Phe.T1.V2.py
# for STS2030+991,forSTS1071
python $rootdir/0_script/Phe.T1.V3.py #need regular
# for NSPT
python $rootdir/0_script/Phe.NSPT.T1.V3.py

```

# step2.phenotype for SF2
```sh
rootdir=/work/home/acd2j8na2s/Work/PRSprofile/
cd ${rootdir}/2_Phenotype
python $rootdir/0_script/SF2.V7.py 'CAS'
```

# step3.GWAS selection for PRS
曙光云计算平台
genome build:hg19
downolad UKBB BMI, SBP, DBP sumstats file
*https://alkesgroup.broadinstitute.org/sumstats_formatted/*
header = SNP,CHR,POS,A1,A2,REF,EAF,Beta,se,P,N,INFO
A1 is EA
A2 is OA
QC:--maf 0.01 --info 0.8
```sh
# GWAS output from UKBB for DBP,SBP,BMI
rootdir=/work/home/acd2j8na2s/Work/PRSprofile/
cd ${rootdir}/3_GWAS
mkdir ${rootdir}/3_GWAS/UKBB
for trait in BMI SBP DBP
do
python $rootdir/0_script/GWAS.UKBB.step1.py \
${rootdir}/3_GWAS/UKBB/${trait}.sumstats \
${rootdir}/3_GWAS/UKBB/${trait}.gwas.QC
done
```
genome build:hg19
download BBJ BMI, SBP, DBP file
header=SNP,CHR,POS,REF,ALT,Frq,Rsq,BETA,SE,P
REF is A2
ALT is A1
QC:--maf 0.01 --Rsq 0.8
BMI: 5644670 SNPs
SBP: 5644670 SNPs
DBP: 5644670 SNPs
```sh
mkdir ${rootdir}/3_GWAS/BBJ
for trait in BMI SBP DBP
do
python $rootdir/0_script/GWAS.UKBB.step1.py \
$rootdir/1_inputdir/GWASoutput/BBJ/${trait}.txt \
$rootdir/1_inputdir/GWASoutput/BBJ/${trait}.gwas.QC
done
```
# step4.merge GWAS files (BBJ or UKBB gwas output) with CAS cohort
retain set of SNPs that overlap between base(gwas output) and target data(all CAS cohort)
## step4.1. 将CAS 中能与EAS匹配的SNP rename 为EASid
前边已完成
## step4.2.只保留CASPMI，STS 和EAS中overloap的SNP(三者SNP数目一样)
前边已完成
## step4.3.进一步保留特定表型-GWAS的SNP，并生成对应的STS,CASPMI,EAS基因型文件
对UKBB-GWAS，1）SNP rename 为EASid；2）保留UKBB-GWAS，CASPMI/STS1071/STS2030 中overloap的SNP(四者SNP数目一样)
对BBJ-GWAS，1）SNP rename 为EASid；2）保留BBJ-GWAS，CASPMI/STS1071/STS2030 中overloap的SNP(四者SNP数目一样)
```sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py36
rootdir=/work/home/acd2j8na2s/Work/PRSprofile
cd ${rootdir}
mkdir $rootdir/4_GWASandData

CASPMIplink=$rootdir/1_Data/CASPMI
STS1071plink=$rootdir/1_Data/STS1071
STS2030plink=$rootdir/1_Data/STS2030
plinkpwd=/work/home/acd2j8na2s/software/plink
for trait in BMI SBP DBP
do
    mkdir $rootdir/4_GWASandData/${trait}
    for db in BBJ UKBB
    do
      mkdir $rootdir/4_GWASandData/${trait}/${db}
      # ---input file
      #/work/home/acd2j8na2s/Work/PRSprofile/3_GWAS/BBJ/BMI.gwas.QC
      db_trait_GWASfile=$rootdir/3_GWAS/${db}/${trait}.gwas.QC
      # ---output file
      gwasrenamefile=$rootdir/4_GWASandData/${trait}/${db}/${trait}.gwas.QC.CAS
      # python to rename GWASid as targetid
      python $rootdir/0_script/GWAS.step4.py ${db_trait_GWASfile} ${CASPMIplink} ${gwasrenamefile}
    done
done

for trait in BMI SBP DBP
do
    for db in BBJ UKBB
    do
    gwasrenamefile=$rootdir/4_GWASandData/${trait}/${db}/${trait}.gwas.QC.CAS
    snptxt=$rootdir/4_GWASandData/${trait}/${db}/${trait}.snp
    awk 'NR>1{print $1}' ${gwasrenamefile} | uniq > $snptxt
    $plinkpwd --bfile ${CASPMIplink} --extract $snptxt --make-bed -out $rootdir/4_GWASandData/${trait}/${db}/CASPMI
    $plinkpwd --bfile ${STS1071plink} --extract $snptxt --make-bed -out $rootdir/4_GWASandData/${trait}/${db}/STS1071
    $plinkpwd --bfile ${STS2030plink} --extract $snptxt --make-bed -out $rootdir/4_GWASandData/${trait}/${db}/STS2030
    done
done  
```

# step5.GWAS meta-analysis
对上一步获得的UKBB-CAS-GWAS 和 BBJ-CAS-GWAS做了meta获得 trans-CASPMI/STS-GWAS
FileNotFoundError: [Errno 2] No such file or directory: '/work/home/acd2j8na2s/Work/PRSprofile/4_GWASandData/BMI/meta/METAANALYSIS1.TBL'

```sh
rootdir=/work/home/acd2j8na2s/Work/PRSprofile
cd ${rootdir}
#wget http://csg.sph.umich.edu/abecasis/metal/download/Linux-metal.tar.gz
#tar -zxvf Linux-metal.tar.gz
metapwd=$rootdir/generic-metal/metal
#metal.txt need to load in windows
metatxt=$rootdir/4_GWASandData/metal.txt

for trait in BMI SBP DBP
do
    rm -r $rootdir/4_GWASandData/${trait}/meta
    mkdir $rootdir/4_GWASandData/${trait}/meta
    cp $rootdir/4_GWASandData/${trait}/BBJ/${trait}.gwas.QC.CAS $rootdir/4_GWASandData/${trait}/meta/BBJ.txt
    cp $rootdir/4_GWASandData/${trait}/UKBB/${trait}.gwas.QC.CAS $rootdir/4_GWASandData/${trait}/meta/UKBB.txt
    cd $rootdir/4_GWASandData/${trait}/meta
    ${metapwd} ${metatxt}
    cd $rootdir/
done

###change metagwas to standform gwas summary file
for trait in BMI SBP DBP
do
metagwas=$rootdir/4_GWASandData/${trait}/meta/METAANALYSIS1.TBL
bbjgwas=$rootdir/4_GWASandData/${trait}/BBJ/${trait}.gwas.QC.CAS
ukbbgwas=$rootdir/4_GWASandData/${trait}/UKBB/${trait}.gwas.QC.CAS
meta2standformgwas=$rootdir/4_GWASandData/${trait}/meta/${trait}.gwas.QC.CAS
python $rootdir/0_script/GWAS.step5.py ${metagwas} ${bbjgwas} ${ukbbgwas} ${meta2standformgwas}
done

#sbatch $rootdir/0_script/test.slurm
#Submitted batch job 3171408
```

# step6.对trans-GWAS，保留trans-GWAS，CASPMI/STS1071/STS2030 中overloap的SNP(四者SNP数目一样)
```sh
for trait in BMI SBP DBP
do
    for db in meta
    do
    gwasrenamefile=$rootdir/4_GWASandData/${trait}/${db}/${trait}.gwas.QC.CAS
    snptxt=$rootdir/4_GWASandData/${trait}/${db}/${trait}.snp
    awk 'NR>1{print $1}' ${gwasrenamefile} | uniq > $snptxt
    $plinkpwd --bfile ${CASPMIplink} --extract $snptxt --make-bed -out $rootdir/4_GWASandData/${trait}/${db}/CASPMI
    $plinkpwd --bfile ${STS1071plink} --extract $snptxt --make-bed -out $rootdir/4_GWASandData/${trait}/${db}/STS1071
    $plinkpwd --bfile ${STS2030plink} --extract $snptxt --make-bed -out $rootdir/4_GWASandData/${trait}/${db}/STS2030
    done
done 
```


# step7.merge all GWAS files (BBJ and UKBB gwas output) with CAS cohort
have rename GWASid as targetid in step4, so here wo want intersection of both
```sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py36
rootdir=/work/home/acd2j8na2s/Work/PRSprofile
cd ${rootdir}
mkdir $rootdir/5_allGWASandData

CASPMIplink=$rootdir/1_Data/CASPMI
STS1071plink=$rootdir/1_Data/STS1071
STS2030plink=$rootdir/1_Data/STS2030
plinkpwd=/work/home/acd2j8na2s/software/plink
for trait in BMI SBP DBP
do
    mkdir $rootdir/5_allGWASandData/${trait}
    mkdir $rootdir/5_allGWASandData/${trait}/plink
    mkdir $rootdir/5_allGWASandData/${trait}/BBJ
    mkdir $rootdir/5_allGWASandData/${trait}/UKBB
    mkdir $rootdir/5_allGWASandData/${trait}/meta
    # ---input file
    gwasrenamefile1=$rootdir/4_GWASandData/${trait}/BBJ/${trait}.gwas.QC.CAS
    gwasrenamefile2=$rootdir/4_GWASandData/${trait}/UKBB/${trait}.gwas.QC.CAS
    gwasrenamefile3=$rootdir/4_GWASandData/${trait}/meta/${trait}.gwas.QC.CAS
    # ---output file
    gwasrenamefile4=$rootdir/5_allGWASandData/${trait}/BBJ/${trait}.gwas.QC.CAS
    gwasrenamefile5=$rootdir/5_allGWASandData/${trait}/UKBB/${trait}.gwas.QC.CAS
    gwasrenamefile6=$rootdir/5_allGWASandData/${trait}/meta/${trait}.gwas.QC.CAS
    # python to rename GWASid as targetid
    python $rootdir/0_script/GWAS.step7.py ${gwasrenamefile1} ${gwasrenamefile2} ${gwasrenamefile3} ${gwasrenamefile4} ${gwasrenamefile5} ${gwasrenamefile6}
done

for trait in BMI SBP DBP
do
    gwasrenamefile=$rootdir/5_allGWASandData/${trait}/BBJ/${trait}.gwas.QC.CAS
    snptxt=$rootdir/5_allGWASandData/${trait}/${trait}.snp
    awk 'NR>1{print $1}' ${gwasrenamefile} | uniq > $snptxt
    $plinkpwd --bfile ${CASPMIplink} --extract $snptxt --make-bed -out $rootdir/5_allGWASandData/${trait}/plink/CASPMI
    $plinkpwd --bfile ${STS1071plink} --extract $snptxt --make-bed -out $rootdir/5_allGWASandData/${trait}/plink/STS1071
    $plinkpwd --bfile ${STS2030plink} --extract $snptxt --make-bed -out $rootdir/5_allGWASandData/${trait}/plink/STS2030
done  
```
#3172702

# step8.GWAS.count.SNP for ST2
```sh
for trait in BMI DBP SBP
do
  for db in UKBB BBJ
  do
  	printf "Trait：",$trait,"db",$db,"\n"
  	if [ ${db} == UKBB ]; then
  		wc -l $rootdir/3_GWAS/UKBB/${trait}.sumstats
  	else
  		wc -l $rootdir/3_GWAS/BBJ/${trait}.txt
  	fi
  	wc -l $rootdir/3_GWAS/$db/${trait}.gwas.QC
  	wc -l $rootdir/4_GWASandData/${trait}/$db/${trait}.gwas.QC.CAS
  	wc -l $rootdir/5_allGWASandData/${trait}/$db/${trait}.gwas.QC.CAS
  done
done

for trait in BMI DBP SBP
do
  for db in meta
  do
  	printf "Trait：",$trait,"db",$db,"\n"
  	wc -l /work/home/acd2j8na2s/Work/PRSprofile/4_GWASandData/${trait}/$db/METAANALYSIS1.TBL
  	wc -l $rootdir/4_GWASandData/${trait}/$db/${trait}.gwas.QC.CAS
  	wc -l $rootdir/5_allGWASandData/${trait}/$db/${trait}.gwas.QC.CAS
  done
done
```

# step9.GWAS plot
曼哈顿图:BBJ UKBB meta
```sh
for trait in BMI DBP SBP
do
	/work/home/acd2j8na2s/soft/R/4.0.3/bin/Rscript $rootdir/0_script/GWAS.cmplot.step9.r $trait
done 
```
曼哈顿图:BBJ UKBB
```sh
for trait in BMI DBP SBP
do
    /work/home/acd2j8na2s/soft/R/4.0.3/bin/Rscript $rootdir/0_script/GWAS.cmplot.step9.V2.r $trait
done 
```

Pmap             输入数据文件
 col             设置不同染色体中点的颜色
 cex             设置点的大小
 pch             设置点的形状
 band            设置不同染色体之间的间隔
 ylim            设置y轴的范围
 bin.size        设置SNP密度图中的窗口大小
 bin.range       设置SNP密度图中图例的范围min，max
 chr.den.col     设置SNP密度的颜色
 cex.axis        设置坐标轴字体的大小
 cex.label        设置坐标轴label字体大小
 lwd.axis        设置坐标轴线的宽度
 plot.type       设置不同的绘图类型，可以设定为 "d", "c", "m", "q" or "b"
 d是snp密度图，c是环形曼哈顿图，m是普通曼哈顿图
 multracks       设置是否需要绘制多个track
 mar             设置图周围白色间隙的大小，应提供4个值，表示底部，左侧，上，右的方向。
 box             是否在曼哈顿图周围加框
 xlab            设置x轴标签
 xticks.pos      设置x刻度标签和x轴之间的距离。
 ylab            设置y轴标签
 ylab.pos        设置y轴label距离轴的距离
 outward         设置点的朝向是否向外
 threshold       设置阈值并添加阈值线，可以为多个，如threshold=c(1,2)
 threshold.col   设置阈值线的颜色,可以设置多个，如threshold.col=c('red','black')
 threshold.lwd   设置阈值线的宽度,如 threshold.lwd=c(1,2)
 threshold.lty   设置阈值线的类型,如 threshold.lty=c(1,2)
 amplify         设置是否放大显著的点
 signal.cex      设置显著点的大小
 signal.pch      设置显著点的形状,可以设置多个,如signal.pch=c(19,19)
 signal.col      设置显著点的颜色
 highlight       设置高光点，如highlight='snp123'
 highlight.cex   设置高光点的大小
 highlight.pch   设置高光点的形状
 highlight.type  设置高光点的类型
 highlight.col   设置高光点的颜色
 highlight.text  设置高光点的文本
 highlight.text.col 设置高光点的文本的颜色
 highlight.text.cex 设置高光点的文本的大小
 highlight.text.xadj设置高光点的文本的水平位置，-1（左）、0（中心）、1（右）
 highlight.text.yadj设置高光点的文本的垂直位置，-1（下）、0（中心）、1（上）
 highlight.text.font设置高光点的文本的字体
 
 
 
 环状曼哈顿
 r               设置圈的半径大小
 chr.labels      设置染色体的标签(密度图和圆曼哈顿图)
 chr.labels.angle设置染色体的标签的角度(密度图和圆曼哈顿图)
 chr.den.col     设置SNP密度图的颜色
 cir.band        设置环状曼哈度图中不同染色体之间的间隔
 cir.chr         设置是否显示染色体的边界
 cir.chr.h       设置染色体边界的高度
 cir.legend      设置是否显示图例
 cir.legend.cex  设置图例字体的大小
 cir.legend.col  设置图例的颜色
 H               设置每个圈的高度
 
 LOG10           设置是否对p-value取log10对数
 conf.int.col    设置QQ图中置信区间的颜色
 file.output     设置是否输出图片
 file            设置输出图片的格式，可以设定为"jpg", "pdf", "tiff"
 dpi             设置输出图片的分辨度
 memo            设置输出图片文件的名字
 height          高
 width           宽
 file            输出文件类型
 file.output     是否输出文件
 main	         标题
 main.cex	     标题大小
 main.font	     标题字体


# step10.GWAS 频率和effect比较图
```sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py36
rootdir=/work/home/acd2j8na2s/Work/PRSprofile
plinkpwd=/work/home/acd2j8na2s/software/plink
all_1kgplink=/work/home/acd2j8na2s/1000G/1000G
EASGfile=/work/home/acd2j8na2s/1000G/EAS_1kg
EURGfile=/work/home/acd2j8na2s/1000G/EUR_1kg
EASsamplefile=/work/home/acd2j8na2s/1000G/EAS.sample
EURsamplefile=/work/home/acd2j8na2s/1000G/EUR.sample
#$plinkpwd --bfile ${all_1kgplink} --keep ${EASsamplefile} --make-bed -out ${EASGfile}
#$plinkpwd --bfile ${all_1kgplink} --keep ${EURsamplefile} --make-bed -out ${EURGfile}
mkdir ${rootdir}/5_allGWASandData/freq
#snptxt1=${rootdir}/5_allGWASandData/BMI/BMI.snp
#snptxt2=${rootdir}/5_allGWASandData/DBP/DBP.snp
#snptxt3=${rootdir}/5_allGWASandData/SBP/SBP.snp
#cat $snptxt1 $snptxt2 $snptxt3 | sort | uniq > ${rootdir}/5_allGWASandData/freq/allsnp
#$plinkpwd --bfile ${EASGfile} --extract ${rootdir}/5_allGWASandData/freq/allsnp --freq -out ${rootdir}/5_allGWASandData/freq/allsnp.EAS.freq 
#$plinkpwd --bfile ${EURGfile} --extract ${rootdir}/5_allGWASandData/freq/allsnp --freq -out ${rootdir}/5_allGWASandData/freq/allsnp.EUR.freq
for trait in BMI DBP SBP
do
    python $rootdir/0_script/GWAS.freqbeta.step10.py $trait
    python $rootdir/0_script/GWAS.freqbeta.step10plot.py $trait
done

for trait in BMI DBP SBP
do
    python $rootdir/0_script/GWAS.freqbeta.step10.V2.py $trait
    python $rootdir/0_script/GWAS.freqbeta.step10plot.V2.py $trait
done

```
# step11.exPRS
```sh
mkdir ${rootdir}/6_ExPRS
mkdir ${rootdir}/6_ExPRS/1.download
cd ${rootdir}/6_ExPRS/1.download
mkdir ${rootdir}/6_ExPRS/2.rePRS
#1. 下载28种ExPRS*5种方法，实际可用19ExPRS
path=${rootdir}/6_ExPRS/1.download
files=$(ls $path)
rm -r ${rootdir}/6_ExPRS/2.rePRS/db.txt
for filename in $files
do
   echo $filename
   awk -v FS="\t" -v OFS="\t" '{if (NR > 13) {print $1,$2,$(NF-2),$(NF-1),$(NF),"'$filename'"}}' ${rootdir}/6_ExPRS/1.download/$filename/ExPRSweb.txt >> ${rootdir}/6_ExPRS/2.rePRS/db.txt
done

head ${rootdir}/6_ExPRS/2.rePRS/db.txt
cat ${rootdir}/6_ExPRS/2.rePRS/db.txt |awk '!a[$6]++{print}'
#2. PRS computation，换一下，用plink --score
#gwasfilr ['CHR', 'SNP', 'BP', 'A1', 'A2', 'BETA']
python $rootdir/0_script/ExPRSandExpouse.step11.1.py 'CASPMI'
python $rootdir/0_script/ExPRSandExpouse.step11.1.py 'STS1071'
python $rootdir/0_script/ExPRSandExpouse.step11.1.py 'STS2030'
#sbatch $rootdir/0_script/test.slurm
#Submitted batch job 3358425,需要看一下结果出来没有
####################error bug#######################################################
find $rootdir/6_ExPRS/2.rePRS/CASPMI/scoredir -name '*.sscore' | wc -l 
#14
find $rootdir/6_ExPRS/2.rePRS/CASPMI/scoredir -name '*.vars' | wc -l 
#16
#LDL Error: --score variant ID 'rs56333512' appears multiple times in main dataset.
#TG Error: --score variant ID 'rs56333512' appears multiple times in main dataset.
grep -w 'rs56333512' $rootdir/1_Data/CASPMI.bim
#7       rs56333512      0       125920259       G       T (not right,del)
#10      rs56333512      0       67715532        A       G
find $rootdir/6_ExPRS/2.rePRS/CASPMI/indexgwasdir -name '*.gwas' | wc -l 
#16
plinkpwd2='/work/home/acd2j8na2s/software/plink2'
#sed 's/要搜索的字符串或正则表达式/替换值/g' 要执行操作的文件名
sed -i 's/rs56333512/rs56333512.del/g' $rootdir/6_ExPRS/2.rePRS/CASPMI/indexgwasdir/LDL.bim
grep -w 'rs56333512' $rootdir/6_ExPRS/2.rePRS/CASPMI/indexgwasdir/LDL.bim
$plinkpwd2 --bfile $rootdir/6_ExPRS/2.rePRS/CASPMI/indexgwasdir/LDL --score $rootdir/6_ExPRS/2.rePRS/CASPMI/indexgwasdir/LDL.gwas 2 4 6 header list-variants cols=scoresums -out $rootdir/6_ExPRS/2.rePRS/CASPMI/scoredir/LDL
sed -i 's/rs56333512/rs56333512.del/g' $rootdir/6_ExPRS/2.rePRS/CASPMI/indexgwasdir/TG.bim
grep -w 'rs56333512' $rootdir/6_ExPRS/2.rePRS/CASPMI/indexgwasdir/TG.bim
$plinkpwd2 --bfile $rootdir/6_ExPRS/2.rePRS/CASPMI/indexgwasdir/TG --score $rootdir/6_ExPRS/2.rePRS/CASPMI/indexgwasdir/TG.gwas 2 4 6 header list-variants cols=scoresums -out $rootdir/6_ExPRS/2.rePRS/CASPMI/scoredir/TG
####################error bug#######################################################
# count snps number
for i in {1..16}
do
    trait=`echo $rootdir/6_ExPRS/filename.txt | awk 'NR=='$i'{print}' $rootdir/6_ExPRS/filename.txt`
    echo $trait
    for pop in CASPMI
    do
        echo $pop
        nodenum=$(cat $rootdir/6_ExPRS/2.rePRS/$pop/scoredir/$trait.sscore.vars | wc -l)
        echo $nodenum
    done
done
#3. 在我们的数据里分析11个ExPRS与expouse的相关性,overloapSNP，PRS in our data,相关性r,p+MSE
python $rootdir/0_script/ExPRSandExpouse.step11.2.py
#4. 在我们的数据里分析可能的暴露因素对应的ExPRS与目标表型的相关性
python $rootdir/0_script/ExPRSandExpouse.step11.3.py
```
###################################### PRS building #########################################
# step12.单一祖先GWAS：C+T，LDpred2，PRSCS，lassosum
# step12.1. provid ld ref for each population
```sh
cd $rootdir
mkdir $rootdir/7_single_ancestry_PRS/
Rpwd=/work/home/acd2j8na2s/soft/R/4.0.3/bin/Rscript

mkdir $rootdir/7_single_ancestry_PRS/LDpred/LDmatrices
cd $rootdir/7_single_ancestry_PRS/LDpred/LDmatrices
wget https://figshare.com/ndownloader/files/37802721

cd $rootdir/
unzip /work/home/acd2j8na2s/Work/PRSprofile/1000-genomes-genetic-maps-master.zip
for i in {1..22}
do
    gunzip /work/home/acd2j8na2s/Work/PRSprofile/1000-genomes-genetic-maps-master/interpolated_OMNI/chr$i.OMNI.interpolated_genetic_map.gz 
done

rootdir=/work/home/acd2j8na2s/Work/PRSprofile
plinkpwd=/work/home/acd2j8na2s/software/plink
all_1kgplink=/work/home/acd2j8na2s/1000G/1000G
EASGfile=/work/home/acd2j8na2s/1000G/EAS_1kg
EURGfile=/work/home/acd2j8na2s/1000G/EUR_1kg
EASsamplefile=/work/home/acd2j8na2s/1000G/EAS.sample
EURsamplefile=/work/home/acd2j8na2s/1000G/EUR.sample
$plinkpwd --bfile ${all_1kgplink} --keep ${EASsamplefile} --make-bed -out ${EASGfile}
$plinkpwd --bfile ${all_1kgplink} --keep ${EURsamplefile} --make-bed -out ${EURGfile}
snptxt1=${rootdir}/5_allGWASandData/BMI/BMI.snp
snptxt2=${rootdir}/5_allGWASandData/DBP/DBP.snp
snptxt3=${rootdir}/5_allGWASandData/SBP/SBP.snp
cat $snptxt1 $snptxt2 $snptxt3 | sort | uniq > ${rootdir}/5_allGWASandData/freq/allsnp
$plinkpwd --bfile ${EASGfile} --extract ${rootdir}/5_allGWASandData/freq/allsnp --make-bed -out ${rootdir}/5_allGWASandData/allsnp.EAS 
$plinkpwd --bfile ${EURGfile} --extract ${rootdir}/5_allGWASandData/freq/allsnp --make-bed -out ${rootdir}/5_allGWASandData/allsnp.EUR
$plinkpwd --bfile ${all_1kgplink} --extract ${rootdir}/5_allGWASandData/freq/allsnp --make-bed -out ${rootdir}/5_allGWASandData/allsnp.1kg

for trait in BMI DBP SBP
do
    mkdir $rootdir/7_single_ancestry_PRS/$trait
    mkdir $rootdir/7_single_ancestry_PRS/$trait/LDmatrices
    mkdir $rootdir/7_single_ancestry_PRS/$trait/LDmatrices/EAS
    mkdir $rootdir/7_single_ancestry_PRS/$trait/LDmatrices/EUR
done

Rpwd=/work/home/acd2j8na2s/soft/R/4.0.3/bin/Rscript
for trait in BMI DBP SBP
do
    $Rpwd $rootdir/0_script/provided-LD.V3.R 'EAS' 'BBJ' $trait
    $Rpwd $rootdir/0_script/provided-LD.V3.R 'EUR' 'UKBB' $trait
done
```
#sbatch $rootdir/0_script/provided-LD.V3.slurm
#Submitted batch job 3839084

# step12.1. ldpred2+lassosum
```sh
#LDpred2
for trait in BMI SBP DBP
do
    mkdir $rootdir/7_single_ancestry_PRS/$trait/ldpred2
    for db in BBJ UKBB meta
    do
        mkdir $rootdir/7_single_ancestry_PRS/$trait/ldpred2/$db
    done
done

plinkpwd=/work/home/acd2j8na2s/software/plink
for trait in BMI SBP DBP
do
  $plinkpwd --bfile $rootdir/5_allGWASandData/$trait/plink/STS2030 \
  --bmerge $rootdir/5_allGWASandData/$trait/plink/STS1071.bed \
  $rootdir/5_allGWASandData/$trait/plink/STS1071.bim \
  $rootdir/5_allGWASandData/$trait/plink/STS1071.fam \
  --make-bed -out $rootdir/5_allGWASandData/$trait/plink/ALLSTS
  if [ ! -f "$rootdir/5_allGWASandData/$trait/plink/ALLSTS-merge.missnp" ];then
  echo "ok1"
  else

  $plinkpwd --bfile $rootdir/5_allGWASandData/$trait/plink/STS2030 \
  --exclude $rootdir/5_allGWASandData/$trait/plink/ALLSTS-merge.missnp \
  --make-bed -out $rootdir/5_allGWASandData/$trait/plink/STS2030.delmissnp

  $plinkpwd --bfile $rootdir/5_allGWASandData/$trait/plink/STS1071 \
  --exclude $rootdir/5_allGWASandData/$trait/plink/ALLSTS-merge.missnp \
  --make-bed -out $rootdir/5_allGWASandData/$trait/plink/STS1071.delmissnp

  $plinkpwd --bfile $rootdir/5_allGWASandData/$trait/plink/STS2030.delmissnp \
  --bmerge $rootdir/5_allGWASandData/$trait/plink/STS1071.delmissnp.bed \
  $rootdir/5_allGWASandData/$trait/plink/STS1071.delmissnp.bim \
  $rootdir/5_allGWASandData/$trait/plink/STS1071.delmissnp.fam \
  --make-bed -out $rootdir/5_allGWASandData/$trait/plink/ALLSTS
  fi
done

for trait in BMI SBP DBP
do
  $plinkpwd --bfile $rootdir/5_allGWASandData/$trait/plink/ALLSTS \
  --bmerge $rootdir/5_allGWASandData/$trait/plink/CASPMI.bed $rootdir/5_allGWASandData/$trait/plink/CASPMI.bim $rootdir/5_allGWASandData/$trait/plink/CASPMI.fam \
  --make-bed -out $rootdir/5_allGWASandData/$trait/plink/ALL
  if [ ! -f "$rootdir/5_allGWASandData/$trait/plink/ALL-merge.missnp" ];then
  echo "ok1"
  else

  $plinkpwd --bfile $rootdir/5_allGWASandData/$trait/plink/ALLSTS \
  --exclude $rootdir/5_allGWASandData/$trait/plink/ALL-merge.missnp \
  --make-bed -out $rootdir/5_allGWASandData/$trait/plink/ALLSTS.delmissnp

  $plinkpwd --bfile $rootdir/5_allGWASandData/$trait/plink/CASPMI \
  --exclude $rootdir/5_allGWASandData/$trait/plink/ALL-merge.missnp \
  --make-bed -out $rootdir/5_allGWASandData/$trait/plink/CASPMI.delmissnp

  $plinkpwd --bfile $rootdir/5_allGWASandData/$trait/plink/ALLSTS.delmissnp \
  --bmerge $rootdir/5_allGWASandData/$trait/plink/CASPMI.delmissnp.bed \
  $rootdir/5_allGWASandData/$trait/plink/CASPMI.delmissnp.bim \
  $rootdir/5_allGWASandData/$trait/plink/CASPMI.delmissnp.fam \
  --pca 6 --make-bed -out $rootdir/5_allGWASandData/$trait/plink/ALL
  fi
done

#sbatch $rootdir/0_script/plinkmerge.slurm
#3584195
```
#产生所有的rds文件
```R 
library(bigsnpr)
library(bigreadr)
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(boot)
on.exit(file.remove('/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/BMI/plink/ALL.rds'), add = TRUE)
on.exit(file.remove('/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/DBP/plink/ALL.rds'), add = TRUE)
on.exit(file.remove('/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/SBP/plink/ALL.rds'), add = TRUE)
on.exit(file.remove('/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/BMI/plink/ALL.bk'), add = TRUE)
on.exit(file.remove('/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/DBP/plink/ALL.bk'), add = TRUE)
on.exit(file.remove('/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/SBP/plink/ALL.bk'), add = TRUE)

snp_readBed('/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/BMI/plink/ALL.bed')
snp_readBed('/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/DBP/plink/ALL.bed')
snp_readBed('/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/SBP/plink/ALL.bed')
```
Error: Not enough disk space to create '/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/SBP/plink/ALL.bk'
df -h ~
#Filesystem      Size  Used Avail Use% Mounted on
#NSCC-XA_work    450G  441G  9.9G  98% /work

#conduct /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/ALL.phe.txt

## 正式的LDpred2+lassosum
```sh
for trait in BMI SBP DBP
do
    echo $trait
    for db in BBJ
    do
        echo $db
        if [ ${db} == BBJ ]; then
            pop=EAS
        else
            pop=EUR
        fi
        echo $pop
        $Rpwd $rootdir/0_script/LDpred2.corr_V2.R $pop $db $trait
  done
done
for trait in BMI SBP DBP
do
    echo $trait
    for db in UKBB
    do
        echo $db
        if [ ${db} == BBJ ]; then
            pop=EAS
        else
            pop=EUR
        fi
        echo $pop
        $Rpwd $rootdir/0_script/LDpred2.corr_V3.R $pop $db $trait
  done
done


for trait in BMI SBP DBP
do
    echo $trait
    for db in BBJ
    do
        echo $db
        if [ ${db} == BBJ ]; then
            pop=EAS
        else
            pop=EUR
        fi
        echo $pop
        $Rpwd $rootdir/0_script/LDpred2.V7.R $rootdir $db $trait $pop
  done
done
#corr用LDpred2.corr.V2.R for each trait in BBJ,  job is running, jobid is 4614544
#sbatch $rootdir/0_script/ldpred2.V7.slurm
#Submitted batch job 4614544

#but in UKBB makes mang 0, and NA or NaN values in the resulting correlation matrix in chr10,6,3, so for UKBB, i try to using download corr, but before computing,  it is need to make sure snps are same in corr and maptest
for trait in BMI SBP DBP
do
    echo $trait
    for db in UKBB
    do
        echo $db
        if [ ${db} == BBJ ]; then
            pop=EAS
        else
            pop=EUR
        fi
        echo $pop
        $Rpwd $rootdir/0_script/LDpred2.V8.R $rootdir $db $trait $pop
  done
done 
#sbatch $rootdir/0_script/ldpred2.V8.slurm
#Submitted batch job 4622139
```
```sh
for trait in BMI SBP DBP
do
    echo $trait
    for db in BBJ UKBB
    do
        echo $db
        if [ ${db} == BBJ ]; then
            pop=EAS
        else
            pop=EUR
        fi
        echo $pop
        $Rpwd $rootdir/0_script/bestLdpred2result.R $rootdir $db $trait $pop
  done
done
```


```sh
for trait in BMI SBP DBP
do
    echo $trait
    for db in BBJ
    do
        echo $db
        if [ ${db} == BBJ ]; then
            pop=EAS
        else
            pop=EUR
        fi
        echo $pop
        $Rpwd $rootdir/0_script/lassosum.V2.R $rootdir $db $trait $pop
  done
done
for trait in BMI SBP DBP
do
    echo $trait
    for db in UKBB
    do
        echo $db
        if [ ${db} == BBJ ]; then
            pop=EAS
        else
            pop=EUR
        fi
        echo $pop
        $Rpwd $rootdir/0_script/lassosum.V3.R $rootdir $db $trait $pop
  done
done
#BMI BBJ lassosum
#sbatch $rootdir/0_script/lassosum.slurm
#Submitted batch job 4621895 ok

#BMI UKBB lassosum
#sbatch $rootdir/0_script/lassosum.slurm
#Submitted batch job 4677106 ？？
address 0x2b7e702f3300, cause 'non-existent physical address'
An irrecoverable exception occurred. R is aborting now ...

#other lassosum for BBJ
#sbatch $rootdir/0_script/lassosumSBPDBPBBJ.slurm
#Submitted batch job 4678158 ??
address 0x2b92a90d1300, cause 'non-existent physical address'
An irrecoverable exception occurred. R is aborting now ...

#other lassosum for UKBB
#sbatch $rootdir/0_script/lassosumSBPDBPUKBB.slurm
#Submitted batch job 4678167 waiting
虽然报错了但是结果输出了？
```
```sh
#进入虚拟环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py36

#需要执行的命令
rootdir=/work/home/acd2j8na2s/Work/PRSprofile
cd ${rootdir}
plinkpwd=/work/home/acd2j8na2s/software/plink
Rpwd=/work/home/acd2j8na2s/soft/R/4.0.3/bin/Rscript

for trait in BMI SBP DBP
do
    echo $trait
    for db in BBJ UKBB
    do
        echo $db
        if [ ${db} == BBJ ]; then
            pop=EAS
        else
            pop=EUR
        fi
        echo $pop
        $Rpwd $rootdir/0_script/bestlassosumresult.R $rootdir $db $trait $pop
  done
done
```

## C+T and SCT
```sh
#test SCT.R
for trait in BMI SBP DBP
do
    mkdir $rootdir/7_single_ancestry_PRS/$trait/sct
    for db in BBJ UKBB
    do
        mkdir $rootdir/7_single_ancestry_PRS/$trait/sct/$db
    done
done
for trait in BMI SBP DBP
do
    echo $trait
    for db in BBJ UKBB
    do
        echo $db
        if [ ${db} == BBJ ]; then
            pop=EAS
        else
            pop=EUR
        fi
        echo $pop
        $Rpwd $rootdir/0_script/SCT.V2.R $rootdir $db $trait $pop
  done
done
#sbatch $rootdir/0_script/SCT.slurm
#Submitted batch job 3910262
```
for BMI UKBB
Error in seq.default(from = log(from), to = log(to), length.out = length.out) : 
  'to' must be a finite number
Calls: snp_grid_PRS ... as.list.environment -> seq_log -> seq -> seq.default
Execution halted


#bestSCTresult.R
```sh
for trait in BMI DBP SBP
do
    echo $trait
    for db in BBJ UKBB
    do
        echo $db
        if [ ${db} == BBJ ]; then
            pop=EAS
        else
            pop=EUR
        fi
        echo $pop
        $Rpwd $rootdir/0_script/bestSCTresult.R $rootdir $db $trait $pop
  done
done
```



## PRScs
```sh
###Getting Started
#Clone this repository using the following git command:
#git clone https://github.com/getian107/PRScs.git
/work/home/acd2j8na2s/Work/PRSprofile/PRScs/PRScs-master/PRScs.py -h

#Download the LD reference panels and extract files:
#LD reference panels constructed using the 1000 Genomes Project phase 3 samples: 
# wget https://personal.broadinstitute.org/hhuang/public/PRS-CSx/Reference/1KG/ldblk_1kg_eas.tar.gz
# wget https://personal.broadinstitute.org/hhuang/public/PRS-CSx/Reference/1KG/ldblk_1kg_eur.tar.gz

#PRScs requires Python packages scipy (https://www.scipy.org/) and h5py (https://www.h5py.org/) installed.
#pip install h5py

##sumstats
for trait in BMI SBP DBP
do
  for db in BBJ UKBB
  do
  $Rpwd ${rootdir}/0_script/PRScs.gwasfile.R \
  $rootdir/5_allGWASandData/$trait/$db \
  $trait.gwas.QC.CAS
  done
done
#SNP is the rs ID, A1 is the effect allele, A2 is the alternative allele, BETA/OR is the effect/odds ratio of the A1 allele, P is the p-value of the effect.

for trait in BMI SBP DBP
do
    mkdir $rootdir/7_single_ancestry_PRS/$trait/PRScs
    for db in BBJ UKBB
    do
        mkdir $rootdir/7_single_ancestry_PRS/$trait/PRScs/$db
    done
done

tar -zxvf $rootdir/PRScs/ldblk_1kg_eas.tar.gz -C $rootdir/PRScs/
tar -zxvf $rootdir/PRScs/ldblk_1kg_eur.tar.gz -C $rootdir/PRScs/

###Using PRS-CS
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
#The output file contains chromosome, rs ID, base position, A1, A2 and posterior effect size estimate for each SNP. 
# Global shrinkage parameter phi will be learnt from the data using a fully Bayesian approach.
sbatch $rootdir/0_script/PRScs.slurm
#Submitted batch job 3822824
for trait in BMI DBP SBP
do
    for db in BBJ UKBB
    do
        cat $rootdir/7_single_ancestry_PRS/$trait/PRScs/${db}_pst_eff_a1_b0.5_phiauto_chr*.txt > $rootdir/7_single_ancestry_PRS/$trait/PRScs/${db}/pst_eff_a1_b0.5_phiauto.txt
    done
done
# PRScs PRS using plink --score
plinkpwd2='/work/home/acd2j8na2s/software/plink2'
#2 4 6
#SNP A1 BETA
for trait in BMI DBP SBP
do
    for db in BBJ UKBB
    do
        for pop in STS2030 CASPMI STS1071
        do
            $plinkpwd2 --bfile ${rootdir}/5_allGWASandData/$trait/plink/$pop --score $rootdir/7_single_ancestry_PRS/$trait/PRScs/${db}/pst_eff_a1_b0.5_phiauto.txt 2 4 6 list-variants cols=scoresums -out $rootdir/7_single_ancestry_PRS/$trait/PRScs/${db}/$pop
            done
        done
done
```
#bestPRscsresult.R
```sh
sbatch /work/home/acd2j8na2s/Work/PRSprofile/0_script/PRSinR.slurm
#4528769
```

## SBayesR
```sh
# download GCTB 
wget https://cnsgenomics.com/software/gctb/download/gctb_2.05beta_Linux.zip
unzip $rootdir/gctb_2.05beta_Linux.zip
gctbpwd=$rootdir/gctb_2.05beta_Linux/gctb

cd $rootdir/gctb_2.05beta_Linux

# step1 download ldm from https://zenodo.org/record/3350914#.ZFhrJqBByyB: GCTB SBayesR shrunk sparse linkage disequilibrium matrices for HM3 variants generated from "Improved polygenic prediction by Bayesian multiple regression on summary statistics" by Lloyd-Jones, Zeng et al. 2019.
wget https://zenodo.org/record/3350914/files/ukbEURu_hm3_sparse.zip?download=1
cat ukbEURu_hm3_sparse.zip?download=1 > ukbEURu_hm3_sparse.zip
unzip ukbEURu_hm3_sparse.zip
cd $rootdir
#sbatch $rootdir/0_script/download.slurm
#

# step2 GCTB summary statistics input format
for trait in BMI SBP DBP
do
  for db in BBJ UKBB
  do
  $Rpwd ${rootdir}/0_script/GCTB.ma_create.R \
  ${rootdir}/5_allGWASandData/$trait/$db \
  $trait.gwas.QC.CAS
  done
done
#sbatch $rootdir/0_script/GCTB.slurm
#ok


#$rootdir/gctb_2.05beta_Linux/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_sparse_mldm_list.txt修改了路径
####################################################################################EAS-LD
mkdir $rootdir/gctb_2.05beta_Linux/EASin1kg_shrunk_sparse
mkdir $rootdir/gctb_2.05beta_Linux/EASplink
EASplink=/work/home/acd2j8na2s/1000G/EAS_1kg
for i in {1..22}
do 
$plinkpwd --bfile ${EASplink} --chr $i --sheep --recode tab --make-bed --out $rootdir/gctb_2.05beta_Linux/EASplink/chr$i
done

for chr in {1..22}
do
  mkdir $rootdir/gctb_2.05beta_Linux/EASin1kg_shrunk_sparse/chr$chr
  t=5000
  nodenum=$(cat $rootdir/gctb_2.05beta_Linux/EASplink/chr$chr.bim | wc -l)
  ((N=$nodenum/$t))
  ((N=$N+1))
  ${gctbpwd} --bfile $rootdir/gctb_2.05beta_Linux/EASplink/chr$chr --make-shrunk-ldm --snp 1-$t --out $rootdir/gctb_2.05beta_Linux/EASin1kg_shrunk_sparse/chr$chr/EASldm_chr$chr

  echo $rootdir/gctb_2.05beta_Linux/EASin1kg_shrunk_sparse/chr$chr/EASldm_chr$chr.snp1-$t.ldm.shrunk > $rootdir/gctb_2.05beta_Linux/EASin1kg_shrunk_sparse/EASldm_chr$chr.mldmlist

  for((i=1;i<=$N;i++))
  do
  ((j=$i))
  ((newt1=$j*$t))
  ((j=$j+1))
  ((newt2=$j*$t))
  echo $newt1-$newt2
  ${gctbpwd} --bfile $rootdir/gctb_2.05beta_Linux/EASplink/chr$chr --make-shrunk-ldm --snp $newt1-$newt2 --out $rootdir/gctb_2.05beta_Linux/EASin1kg_shrunk_sparse/chr$chr/EASldm_chr$chr

  echo $rootdir/gctb_2.05beta_Linux/EASin1kg_shrunk_sparse/chr$chr/EASldm_chr$chr.snp$newt1-$newt2.ldm.shrunk >> $rootdir/gctb_2.05beta_Linux/EASin1kg_shrunk_sparse/EASldm_chr$chr.mldmlist

  done
done

for chr in {1..22}
do
cat $rootdir/gctb_2.05beta_Linux/EASin1kg_shrunk_sparse/EASldm_chr$chr.mldmlist
# Merge the LD matrix chunks into a single file
${gctbpwd} --mldm $rootdir/gctb_2.05beta_Linux/EASin1kg_shrunk_sparse/EASldm_chr$chr.mldmlist --make-shrunk-ldm --out $rootdir/gctb_2.05beta_Linux/EASin1kg_shrunk_sparse/EASldm_chr$chr
done
#sbatch $rootdir/0_script/GCTB.step2.slurm
#Submitted batch job 

####################################################################################
for trait in BMI SBP DBP
do
    mkdir $rootdir/7_single_ancestry_PRS/$trait/SBayesR
    for db in BBJ UKBB
    do
        mkdir $rootdir/7_single_ancestry_PRS/$trait/SBayesR/$db
    done
done

# step3 Calling SBayesR
for trait in BMI SBP DBP
do
  for db in BBJ
  do
  mkdir $rootdir/2_PRS/$trait/$db/SBayesR
  ${gctbpwd} --sbayes R \
  --mldm $rootdir/gctb_2.05beta_Linux/EASin1kg_shrunk_sparse/EASin1kg_shrunk_sparse.txt \
  --pi 0.95,0.02,0.02,0.01 \
  --gamma 0.0,0.01,0.1,1 \
  --gwas-summary ${rootdir}/5_allGWASandData/$trait/$db/$trait.gwas.QC.CAS.ma \
  --chain-length 10000 \
  --burn-in 2000 \
  --out-freq 10 \
  --out $rootdir/7_single_ancestry_PRS/$trait/SBayesR/$db/SBayesR.output 2>&1 | tee "$rootdir/7_single_ancestry_PRS/$trait/SBayesR/$db/SBayesR.recode.log"
  done
done


for trait in BMI SBP DBP
do
  for db in UKBB
  do
  mkdir $rootdir/2_PRS/$trait/$db/SBayesR
  ${gctbpwd} --sbayes R \
  --mldm $rootdir/gctb_2.05beta_Linux/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_sparse_mldm_list.txt \
  --pi 0.95,0.02,0.02,0.01 \
  --gamma 0.0,0.01,0.1,1 \
  --gwas-summary ${rootdir}/5_allGWASandData/$trait/$db/$trait.gwas.QC.CAS.ma \
  --chain-length 10000 \
  --burn-in 2000 \
  --out-freq 10 \
  --out $rootdir/7_single_ancestry_PRS/$trait/SBayesR/$db/SBayesR.output 2>&1 | tee "$rootdir/7_single_ancestry_PRS/$trait/SBayesR/$db/SBayesR.recode.log"
  done
done
```

## PRS-csx
```sh
for trait in BMI SBP DBP
do
    mkdir $rootdir/7_single_ancestry_PRS/$trait/PRScsx
done
###Using PRS-CS
#ref_dir：LD参考面板的路径，路径下应包含相应群体的参考面板以及snp list. 例如，纳入群体为EUR以及EAS，指定路径为：./ldref ，那么该路径下应该有 ldblk_1kg_eas，ldblk_1kg_eur 这两个文件夹， 以及snpinfo_mult_1kg_hm3这个文件。
#bim_prefix：目标数据集的bim文件。
#sst_file：sumstats的完整路径，由逗号分隔。
#n_gwas：sumstats的样本量大小，由逗号分隔，顺序与SUM_STATS_FILE一致。
#pop：对应的群体，可以为 AFR, AMR, EAS, EUR, SAS，由逗号分隔，顺序与SUM_STATS_FILE一致。
#out_dir: 输出的路径
#META_FLAG ： 如果为True，则输出inverse-variance-weighted meta-analysis of the population-specific posterior effect size estimates。
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

for trait in BMI DBP SBP
do
    cat $rootdir/7_single_ancestry_PRS/$trait/PRScsx/PRScsx_META_pst_eff_a1_b0.5_phiauto_chr*.txt \
    > $rootdir/7_single_ancestry_PRS/$trait/PRScsx/PRScsx_META_pst_eff_a1_b0.5_phiauto.txt
done
# PRScsx PRS using plink --score
plinkpwd2='/work/home/acd2j8na2s/software/plink2'
#2 4 6
#SNP A1 BETA
for trait in BMI DBP SBP
do
    for pop in STS2030 CASPMI STS1071
        do
            mkdir $rootdir/7_single_ancestry_PRS/$trait/PRScsx/score
            $plinkpwd2 --bfile ${rootdir}/5_allGWASandData/$trait/plink/$pop --score $rootdir/7_single_ancestry_PRS/$trait/PRScsx/PRScsx_META_pst_eff_a1_b0.5_phiauto.txt 2 4 6 list-variants cols=scoresums -out $rootdir/7_single_ancestry_PRS/$trait/PRScsx/score/$pop
        done
done
```
#bestPRscsresult.R
```sh
sbatch /work/home/acd2j8na2s/Work/PRSprofile/0_script/PRSinRforPRScsx.slurm
#4529149
```
#sbatch $rootdir/0_script/PRScsx.slurm
#Submitted batch job 3925131
1）UnboundLocalError: local variable 'ref_dict' referenced before assignment
cp /work/home/acd2j8na2s/Work/PRSprofile/PRScs/snpinfo_mult_1kg_hm3.txt /work/home/acd2j8na2s/Work/PRSprofile/PRScs/snpinfo_mult_1kg_hm3
OK


## CT-SLEB
*https://andrewhaoyu.github.io/CTSLEB/*
```sh
for trait in BMI SBP DBP
do
    mkdir $rootdir/7_single_ancestry_PRS/$trait/CTSLEB
done

for trait in BMI SBP DBP
do
    $Rpwd $rootdir/0_script/CT-SLEB.V4.R $trait
done

sbatch $rootdir/0_script/CT-SLEB.slurm
#4246842
```

## jointPRS
*JointPRS, https://github.com/LeqiXu/JointPRS
```sh
source ~/miniconda3/etc/profile.d/conda.sh

#conda activate py36
#conda deactivate
#1. JointPRS Installation
cd /work/home/acd2j8na2s/software
git clone https://github.com/LeqiXu/JointPRS.git
cd /work/home/acd2j8na2s/software/JointPRS
conda env create -f environment.yml
conda activate JointPRS
python setup.py build_ext --inplace

#2. LD Reference Panel Download
# use reference panels from PRS-CSx and you can follow their instructions to download them. It is strongly recommended to create two subfolders within your reference directory:
#have download in /work/home/acd2j8na2s/Work/PRSprofile/PRScs

#3. Summary Statistics Preparation
#header= SNP               A1      A2      BETA            P
#A1: the effect allele.
#In addition, you need to obtain the sample size for the summary statistics, and take the median value if the sample size is different across SNPs.
#${rootdir}/5_allGWASandData/$trait/BBJ/$trait.gwas.QC.CAS.PRScs.txt

#4. JointPRS Model Implementation
sbatch $rootdir/0_script/jointPRS.slurm
#19994329
#squeue
cat ${rootdir}/slurm-19994329.out

sbatch $rootdir/0_script/jointPRS.SBP.slurm
cat ${rootdir}/slurm-20682796.out
```
param_list = c("auto","1e-06","1e-04","1e-02","1e+00")
耗时太久，一个表型需要7-9天时间。还是只是用auto。

```sh
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
    cat $rootdir/7_single_ancestry_PRS/$trait/JointPRS/JointPRS_auto_EAS_EUR_11_1KG_EAS_pst_eff_a1_b0.5_phiauto_chr*.txt \
    > $rootdir/7_single_ancestry_PRS/$trait/JointPRS/JointPRS_auto_EAS_phiauto.txt
done
# jointPRS using R
for trait in BMI DBP SBP
do
    $Rpwd /work/home/acd2j8na2s/Work/PRSprofile/0_script/PRSinRforJointPRS.R $trait
done
```

## ROSPER
*https://github.com/Jingning-Zhang/PROSPER
```sh
#进入虚拟环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py36

#step1.install
cd /work/home/acd2j8na2s/software
git clone https://github.com/wkentaro/gdown.git
cd gdown
pip install gdown

cd /work/home/acd2j8na2s/software
git clone https://github.com/Jingning-Zhang/PROSPER.git
## go to the directory
cd /work/home/acd2j8na2s/software/PROSPER

#step2.download reference using the 1000 Genomes Project phase 3 samples for variants in either Hapmap3 or MEGA.
#EUR reference (~8.6G); 
#EAS reference (~6.2G); 
#ref_bim.txt
## download reference SNPs
gdown 1PtD4qk7EBPxdhkGrKrG8OktKxaSDF9AT #error,download in wins
## download EUR reference LD
gdown 1ger1-jsoD73vCez5vN6h4QD5TMdU_WS6 #error,download in wins
tar -zxvf EUR.tar.gz
## download EAS reference LD
gdown 1NzltrpebQiaHYRIN67nTqmRChJ_0MEze #error,download in wins
tar -zxvf EAS.tar.gz

## You will still need to install plink2 and the required R packages
#Launch R and install required libraries:
#install.packages(c('optparse','bigreadr','readr','stringr', 'caret', 'SuperLearner', 'glmnet', 'MASS', 'Rcpp', 'RcppArmadillo', 'inline', 'doMC', 'foreach'))

#需要执行的命令
rootdir=/work/home/acd2j8na2s/Work/PRSprofile
Rpwd=/work/home/acd2j8na2s/soft/R/4.0.3/bin/Rscript
cd ${rootdir}
for trait in BMI DBP SBP
do
mkdir ${rootdir}/7_single_ancestry_PRS/$trait/ROSPER
done
#x=/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/BMI/LDmatrices

#Using PROSPER
#There are two scripts PROSPER.R and tuning_testing.R. 
#The first one is used to get solutions from all tuning parameters setting, 
#The second one is used to perform tuning/testing and get the final ensembled PRS. 
#The parameters are explained in this section. 
PROSPER.R --PATH_package --PATH_out --FILE_sst --pop --lassosum_param --chrom --Ll --Lc --verbose --NCORES

tuning_testing.R --PATH_out --PATH_plink --prefix --SL_library --linear_score --bfile_tuning --pheno_tuning --covar_tuning --testing --bfile_testing --pheno_testing --covar_testing --verbose --cleanup --NCORES
#PATH_package (required): Full path to the directory mentioned above that contains: 1. a folder scripts downloaded from github 2. LD reference panels downloaded from Google Drive.
#PATH_plink (required): Full path to plink2 executable. Please do not use plink1.9
#PATH_out (required): Full path to the directory for output results
#FILE_sst (required): Full path and the file name of the GWAS summary statistics from multiple populations, separated by comma. An example of the format, rsid           chr     a1     a0    beta       beta_se    n_eff
#a1: effective allele

plinkpwd2='/work/home/acd2j8na2s/software/plink2'

#the summary statistics files are suggested to be cleaned as follows before using:
##The rsid of variants in reference panels can be found in ref_bim.txt
#Keep variants in reference panels
#Alleles are suggested to match to reference panels to avoid flipping strands
#For each population, remove variants whose allele frequencies are less than 1%. If your summary statistics does not contain information of population-specific allele frequencies, we recommanded to use that in 1000G, which can be found in ref_af.txt.
#Remove variants with small sample size (90% of median sample size per variant).
python $rootdir/0_script/PROSPER.clean.py
#need optimal tuning parameters in lassosum2
#必须cd到这里cd /work/home/acd2j8na2s/software/PROSPER
#creates a directory $PATH_out/before_ensemble/
#writes two files score_file.txt and score_param.txt inside it.
#score_file.txt contains all solutions from PROSPER
#score_param.txt contains their ancestry origin and corresponding tuning parameter settings
cd /work/home/acd2j8na2s/software/PROSPER
sbatch $rootdir/0_script/PROSPER.V2.slurm
#20593890
cat /work/home/acd2j8na2s/software/PROSPER/slurm-20593890.out
```
#Error in { : 
#  task 15 failed - "number of columns of matrices must match (see arg 2)"
#Calls: %dopar% -> <Anonymous>
#Execution halted
#没有15号染色体?why
##ERROR: 15号染色体 in ## Step 2.4. Clean PRSs into a matrix (#variant X #grid_search)
/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/BMI/PROSPER/1/before_ensemble/score_file.txt
#rsid    a1      a0      score1-50
#42507 SNPs for chr1
/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/BMI/PROSPER/1/before_ensemble/score_param.txt
#score_origin    delta_EAS       delta_EUR       lambda_EAS      lambda_EUR      c       sparsity_nonzero_percentage

```sh
package='/work/home/acd2j8na2s/software/PROSPER'
path_plink='/work/home/acd2j8na2s/software/plink2'
target_pop='EAS'
Rpwd=/work/home/acd2j8na2s/soft/R/4.0.3/bin/Rscript
for trait in BMI DBP SBP
do
    path_result=/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/$trait/
    $Rpwd ${package}/scripts/tuning_testing.R \
    --PATH_plink ${path_plink} \
    --PATH_out ${path_result}/PROSPER \
    --prefix ${target_pop} \
    --testing TRUE \
    --bfile_tuning /work/home/acd2j8na2s/Work/PRSprofile/4_GWASandData/$trait/meta/STS2030 \
    --pheno_tuning /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/STS2030.$trait.pheno.fam \
    --bfile_testing /work/home/acd2j8na2s/Work/PRSprofile/4_GWASandData/$trait/meta/CASPMI \
    --pheno_testing /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/CASPMI.$trait.pheno.fam \
    --cleanup F \
    --NCORES 5
done

#head /work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/BMI/PROSPER/after_ensemble_EAS/PROSPER_prs_file.txt
#rsid    a1      a0      weight
#478776 SNPs

#compute PRS for STS2030,CASPMI and STS1071 with R2
#1st.beta文件增加chr,BP，'chr','rsid','pos','a1','a0','beta'
python PRSinRforPROSPER.py
for trait in BMI DBP SBP
do
    $Rpwd /work/home/acd2j8na2s/Work/PRSprofile/0_script/PRSinRforPROSPER.R $trait
done
```



## Weighted-PRS
```sh
$Rpwd WeightedPRS.V1.R 

```

## PolyPred
```sh
$rootdir/0_script/ploypred.test.BMI.md
$rootdir/0_script/ploypred.test.SBP.md
$rootdir/0_script/ploypred.test.DBP.md
```

ALL.result.compare.V2.csv
```R
library('ggplot2')
library(tidyr)
library("tidyverse")
#library(ggthemes)
#library(RColorBrewer)
#library(colorspace)

setwd('E:\\04Zhangyaning\\PRS\\2_PRS_V3\\7_single_ancestry_PRS')
mydata <- read.table('ALL.result.compare.V2.csv',sep=',',header=T)
mydata1 <- mydata[c('Trait', 'Method', 'GWAS', 'TrainingR2','TrainingR2_ci005','TrainingR2_ci095')]
mydata1['population']='Tuning'
mydata2 <- mydata[c('Trait', 'Method', 'GWAS', 'ValidationR2','ValidationR2_ci005','ValidationR2_ci095')]
mydata2['population']='Testing'
mydata3 <- mydata[c('Trait', 'Method', 'GWAS', 'TestingR2','TestingR2_ci005','TestingR2_ci095')]
mydata3['population']='Validation'
names(mydata1) <- c('Trait', 'Method', 'GWAS', 'R2','R2_CI5','R2_CI95','population')
names(mydata2) <- c('Trait', 'Method', 'GWAS', 'R2','R2_CI5','R2_CI95','population')
names(mydata3) <- c('Trait', 'Method', 'GWAS', 'R2','R2_CI5','R2_CI95','population')
newz <- rbind(mydata1,mydata2,mydata3)

newz <- newz %>%  mutate(population = factor(population, levels = c("Tuning","Testing","Validation")))

newz <- newz %>%
  mutate(GWAS = factor(GWAS, levels = c("BBJ", "UKBB", "Multi-ancestry")))

newz <- unite(newz,"GWASMethod",c("GWAS","Method"), sep="-", remove = F)

newz <- newz %>%
  mutate(GWASMethod = factor(GWASMethod, levels = c("BBJ-C+T","BBJ-SCT","BBJ-PRScs","BBJ-ldpred2","BBJ-lassosum","UKBB-C+T","UKBB-SCT", "UKBB-PRScs", "UKBB-ldpred2","UKBB-lassosum","Multi-ancestry-PRScsx","Multi-ancestry-CT-SLEB","Multi-ancestry-PolyPredP+","Multi-ancestry-JointPRS","Multi-ancestry-PROSPER")))
levels(newz$GWASMethod)

#xx <- subset(newz,Trait=="BMI" & population=="Validation")

p <- ggplot(newz, aes(x = GWASMethod, y = R2, fill=GWASMethod)) + 
    geom_bar(color = 'black',position=position_dodge(), stat="identity") +
    scale_fill_manual(values = c('#C6E2FF','#87CEEB','#7EC0EE','#6CA6CD','#104E8B','#CDC5BF','#8B8682','#8B7765','#787878','#3D3D3D','#FFC1C1','#FA8072','#FF4500','#FF4040','#8B0000'))+
    #geom_errorbar(aes(ymin=R2_CI5, ymax=R2_CI95), width=.4,    position=position_dodge(.9))+
    theme_bw() +
    theme( 
        #legend.position="none",
        axis.ticks = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.title.x=element_text(face="bold", size=35, color = "black", vjust = -0.5),  
        axis.title.y=element_text(face="bold", size=35, color = "black"),  
        axis.text.x = element_blank(), 
        axis.text.y=element_text(size=15, face="bold", color = "black")) +
    scale_y_continuous(limits = c(0, 0.12),breaks = seq(0,0.12,0.02))

p + facet_grid(Trait ~ population)+
    theme(
      strip.text.x = element_text(
        size = 20, color = "black", face = "bold"
        ), # 这里设置x轴方向的字体类型，
      strip.text.y = element_text(
        size = 20, color = "black", face = "bold"
        ) # 这里设置y轴方向的字体类型，
      ) +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave('E:\\04Zhangyaning\\PRS\\2_PRS_V3\\7_single_ancestry_PRS\\allcohort.bestPRS.tiff', plot = last_plot(),
  scale = 1, width = 50, height = 25, units =c("cm"),
  dpi = 300, limitsize = TRUE)

setEPS()
postscript("E:\\04Zhangyaning\\PRS\\2_PRS_V3\\7_single_ancestry_PRS\\allcohort.bestPRS.eps")
dev.off()

```

```R
library('ggplot2')
library(tidyr)
library("tidyverse")
#library(ggthemes)
#library(RColorBrewer)
#library(colorspace)

setwd('E:\\04Zhangyaning\\PRS\\2_PRS_V3\\7_single_ancestry_PRS')
mydata <- read.table('ALL.result.compare.V2.csv',sep=',',header=T)
mydata1 <- mydata[c('Trait', 'Method', 'GWAS', 'TrainingR2','TrainingR2_ci005','TrainingR2_ci095')]
mydata1['population']='Tuning'
mydata2 <- mydata[c('Trait', 'Method', 'GWAS', 'TestingR2','TestingR2_ci005','TestingR2_ci095')]
mydata2['population']='Testing'
mydata3 <- mydata[c('Trait', 'Method', 'GWAS', 'ValidationR2','ValidationR2_ci005','ValidationR2_ci095')]
mydata3['population']='Validation'
names(mydata1) <- c('Trait', 'Method', 'GWAS', 'R2','R2_CI5','R2_CI95','population')
names(mydata2) <- c('Trait', 'Method', 'GWAS', 'R2','R2_CI5','R2_CI95','population')
names(mydata3) <- c('Trait', 'Method', 'GWAS', 'R2','R2_CI5','R2_CI95','population')
newz <- rbind(mydata2,mydata3)
newz <- newz %>%  mutate(population = factor(population, levels = c("Testing","Validation")))

newz <- newz %>%
  mutate(GWAS = factor(GWAS, levels = c("BBJ", "UKBB", "Multi-ancestry")))

newz <- unite(newz,"GWASMethod",c("GWAS","Method"), sep="-", remove = F)

newz <- newz %>%
  mutate(GWASMethod = factor(GWASMethod, levels = c("BBJ-C+T","BBJ-SCT","BBJ-PRScs","BBJ-ldpred2","BBJ-lassosum","UKBB-C+T","UKBB-SCT", "UKBB-PRScs", "UKBB-ldpred2","UKBB-lassosum","Multi-ancestry-PRScsx","Multi-ancestry-CT-SLEB","Multi-ancestry-PolyPredP+","Multi-ancestry-JointPRS","Multi-ancestry-PROSPER")))
levels(newz$GWASMethod)

#xx <- subset(newz,Trait=="BMI" & population=="Validation")

p <- ggplot(newz, aes(x = GWASMethod, y = R2, fill=GWASMethod)) + 
    geom_bar(color = 'black',position=position_dodge(), stat="identity") +
    scale_fill_manual(values = c('#C6E2FF','#87CEEB','#7EC0EE','#6CA6CD','#104E8B','#CDC5BF','#8B8682','#8B7765','#787878','#3D3D3D','#FFC1C1','#FA8072','#FF4500','#FF4040','#8B0000'))+
    #geom_errorbar(aes(ymin=R2_CI5, ymax=R2_CI95), width=.4,    position=position_dodge(.9))+
    theme_bw() +
    theme( 
        #legend.position="none",
        axis.ticks = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.title.x=element_text(face="bold", size=35, color = "black", vjust = -0.5),  
        axis.title.y=element_text(face="bold", size=35, color = "black"),  
        axis.text.x = element_blank(), 
        axis.text.y=element_text(size=15, face="bold", color = "black")) +
    scale_y_continuous(limits = c(0, 0.1),breaks = seq(0,0.1,0.02))

p + facet_grid(population ~ Trait)+
    theme(
      strip.text.x = element_text(
        size = 20, color = "black", face = "bold"
        ), # 这里设置x轴方向的字体类型，
      strip.text.y = element_text(
        size = 20, color = "black", face = "bold"
        ) # 这里设置y轴方向的字体类型，
      ) +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave('E:\\04Zhangyaning\\PRS\\2_PRS_V3\\7_single_ancestry_PRS\\testing.bestPRS.tiff', plot = last_plot(),
  scale = 1, width = 50, height = 25, units =c("cm"),
  dpi = 300, limitsize = TRUE)

setEPS()
postscript("E:\\04Zhangyaning\\PRS\\2_PRS_V3\\7_single_ancestry_PRS\\testing.bestPRS.eps")
dev.off()
```
# analysis bestPRS in disease prediction,obesity;hypertension;obesity and hypertension
```sh
#需要执行的命令
rootdir=/data1/yaning/work/PRSprofile
cd ${rootdir}
Rpwd=/usr/bin/Rscript

#compute OR/SD and Q2-5 OR: BMIPRS for fat; (DBPPRS+SBPPRS)/2 for HTN;
python /data1/yaning/work/PRSprofile/0_script/multinom.step1.V5.py
$Rpwd /data1/yaning/work/PRSprofile/0_script/multinom.step2.V4.R
```


#AUC+OR+boxplot组合图，V6是正式的5折参考loocv
```sh
#python /data1/yaning/work/PRSprofile/0_script/FigbestPRS.OR.AUC.V7.py STS2030
#python /data1/yaning/work/PRSprofile/0_script/FigbestPRS.OR.AUC.V7.py CASPMI
#python /data1/yaning/work/PRSprofile/0_script/FigbestPRS.OR.AUC.V7.py STS1071
```

#PRS_BMI, PRS_HTN crosstable
```sh
#python /work/home/acd2j8na2s/Work/PRSprofile/0_script/crosstable.py STS1071 未完成，先等下
```

#OrH AUC图
```sh
#python /data1/yaning/work/PRSprofile/0_script/FigbestPRS.OrH.AUC.V3.py STS2030
#python /data1/yaning/work/PRSprofile/0_script/FigbestPRS.OrH.AUC.V3.py CASPMI
#python /data1/yaning/work/PRSprofile/0_script/FigbestPRS.OrH.AUC.V3.py STS1071
```
```sh
#NRIforM0andM1,5%,阈值
#python /work/home/acd2j8na2s/Work/PRSprofile/0_script/bestPRS.allsample.NRI.py STS2030
#python /work/home/acd2j8na2s/Work/PRSprofile/0_script/bestPRS.allsample.NRI.V2.py STS2030
#python /work/home/acd2j8na2s/Work/PRSprofile/0_script/bestPRS.allsample.NRI.V2.py CASPMI
#python /work/home/acd2j8na2s/Work/PRSprofile/0_script/bestPRS.allsample.NRI.V2.py STS1071


#NRI:sex
# del python /work/home/acd2j8na2s/Work/PRSprofile/0_script/bestPRS.allsample.NRI.V2.sexdiff.py STS1071
```

#PRS新的主图，prsboxplot&OR&AUC for fat, HTN and OrH
```sh
$Rpwd /data1/yaning/work/PRSprofile/0_script/FigbestPRS.all.step2.V3.R
python /data1/yaning/work/PRSprofile/0_script/FigbestPRS.all.AUC.py STS1071
python /data1/yaning/work/PRSprofile/0_script/FigbestPRS.all.AUC.py STS2030
python /data1/yaning/work/PRSprofile/0_script/FigbestPRS.all.AUC.py CASPMI
```
#计算r2&95%CI
```sh
Rpwd=/work/home/acd2j8na2s/soft/R/4.0.3/bin/Rscript
$Rpwd $rootdir/0_script/RwithandwithoutsexandageforPRS.R
```
#95%CI AUC
```R
library(boot)
library(pROC)
# boot的使用方式很奇怪
get_auc <- function(data, ind, outcome, predictor){
  d = data[ind,] #这句必须加
  au <- as.numeric(auc(pROC::roc(d[,outcome], d[,predictor],quiet=T)))
  au
}
#HTN_definition
pop<-'STS1071'
for (target in c('OrH')) {
    print(target)
    phenofilefile <- paste('E:/04Zhangyaning/PRS/2_PRS_V3/8_bestPRS/',pop,target,'.prob.csv',sep = "")
    dt = read.csv(phenofilefile, header = T, row.names = 1)
    get_auc(dt, outcome="target",predictor="M1_prob.1")
    set.seed(123)
    ba <- boot(dt, get_auc, R = 1000,
           outcome="target",predictor="M1_prob.1")
    print(ba)
    print(boot.ci(ba,conf = 0.95, type = "perc"))
}
```
#PRS for sex boxplot
```py
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np 
import os 
import seaborn as sns 
from scipy import stats
import statsmodels.api as sm
from  scipy.stats import chi2_contingency
import numpy as np
from scipy.stats import t
import statsmodels.formula.api as smf
import sys
from sklearn.preprocessing import StandardScaler
import matplotlib.font_manager as font_manager
import matplotlib as mpl
import random
import shutil
from functools import reduce
from statsmodels.formula.api import ols
from scipy import interp
import scipy.stats as st
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LinearRegression,LogisticRegression
from sklearn.model_selection import LeaveOneOut,StratifiedKFold,KFold
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.metrics import median_absolute_error, r2_score,confusion_matrix,explained_variance_score,mean_absolute_error
from statannotations.Annotator import Annotator
#pip install statannotations
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import KFold,StratifiedKFold

plt.switch_backend('agg')
plt.cla()
plt.close("all")
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300

def getPRSadjage(data,PRSname):
    X = data['age']
    y = data[PRSname]
    X = sm.add_constant(X)
    model = sm.OLS(y, X)
    results = model.fit()
    residuals = results.resid
    return(residuals)

font_size = 30
plt.rcParams['font.size'] = font_size
plt.rcParams['figure.figsize'] = (12,6)
mpl.rc('axes', lw=5)
workdir='/data1/yaning/work/PRSprofile/8_bestPRS'
fig = plt.figure()
#subplot(numRows, numCols, plotNum)
axes = fig.subplots(1,2)

df=pd.read_csv('/data1/yaning/work/PRSprofile/8_bestPRS/STS1071OrH.prob.csv')
dat1=df[['id','sex','PRS_BMI']]
dat2=df[['id','sex','PRS_HTN']]
dat1['PRS_BMIadj'] = getPRSadjage(df,'PRS_BMI')
dat2['PRS_HTNadj'] = getPRSadjage(df,'PRS_HTN')


ax1 = axes[0]
prsscorename='PRS_BMIadj'
c1, c0, c2 = sns.color_palette('Set1', 3)
dist1 = dat1[dat1['sex']==1][prsscorename]
dist0 = dat1[dat1['sex']==2][prsscorename]
mean1 = dist1.mean()
std1 = dist1.std()
mean0 = dist0.mean()
std0 = dist0.std()
print(mean0,std0)
print(mean1,std1)
stat_val01, p_val01 = stats.ttest_ind(dist0, dist1, equal_var=False)
p_val01 = "{:.2e}".format(p_val01)
print(p_val01)
my_pal = {1: "#f9766e", 2:"lightgreen"}
sns.boxplot(data=dat1, x='sex', y=prsscorename,ax=ax1,linewidth=5,palette=my_pal,width=0.8)
ax1.set_xticklabels(['Male','Female'],rotation = 0,fontsize = font_size)
pairs=[(1,2)]
annotator = Annotator(ax1, pairs, data=dat1, x='sex', y=prsscorename)
annotator.configure(test='Mann-Whitney', text_format='star',line_height=0.03,line_width=5)
annotator.apply_and_annotate()
ax1.set_ylabel('PRS_BMI',fontsize=font_size,fontweight='bold')
ax1.set_xlabel('Sex',fontsize=font_size,fontweight='bold')
ax1.set_ylim([-3,6])

ax2 = axes[1]
prsscorename='PRS_HTNadj'
c1, c0, c2 = sns.color_palette('Set1', 3)
dist1 = dat2[dat2['sex']==1][prsscorename]
dist0 = dat2[dat2['sex']==2][prsscorename]
mean1 = dist1.mean()
std1 = dist1.std()
mean0 = dist0.mean()
std0 = dist0.std()
print(mean0,std0)
print(mean1,std1)
stat_val01, p_val01 = stats.ttest_ind(dist0, dist1, equal_var=False)
p_val01 = "{:.2e}".format(p_val01)
print(p_val01)
my_pal = {1: "#f9766e", 2:"lightgreen"}
sns.boxplot(data=dat2, x='sex', y=prsscorename,ax=ax2,linewidth=5,palette=my_pal,width=0.8)
ax2.set_xticklabels(['Male','Female'],rotation = 0,fontsize = font_size)
pairs=[(1,2)]
annotator = Annotator(ax2, pairs, data=dat2, x='sex', y=prsscorename)
annotator.configure(test='Mann-Whitney', text_format='star',line_height=0.03,line_width=5)
annotator.apply_and_annotate()
ax2.set_ylabel('PRS_HTN',fontsize=font_size,fontweight='bold')
ax2.set_xlabel('Sex',fontsize=font_size,fontweight='bold')
ax2.set_ylim([-3,6])
plt.tight_layout()
plt.savefig('/data1/yaning/work/PRSprofile/8_bestPRS/PRS.sexplot.png',bbox_inches='tight')
```



#MRS for BMI using 187 CpGs, at P<1x10-7 in D, and <0.05 in R and with directional consistency
#MRS for BP using 34 cross-validated CpG sites from the meta-analysis
沙箱可以了。
操作流程是，先网页登录https://box.biosino.org:4430/dashboard进入堡垒机系统，账号:zhangyaning，密码: Picb@2024+，随后通过web或者ssh途径登录沙箱服务器
数据存放在了/home/zhangyaning/data/目录下，因为home目录空间比较小，如果师姐你要存东西或者跑程序生成output，可以在/data/run/目录下运行
```sh
cd /data/run/

#需要执行的命令
rootdir=/data/run/
cd ${rootdir}
#plinkpwd=/work/home/acd2j8na2s/software/plink
Rpwd=/usr/bin/R
#step0. get cpgfile from EWAS Catalog
#python $rootdir/0_script/MRS.step0.getcpgfile.V2.py
$Rpwd $rootdir/0_script/MRS.step1.extractcpg.V3.R

python $rootdir/0_script/MRS.step2.splittraintest.V2.py
#step2.compute beta using 弹性网:
python $rootdir/0_script/MRS.step2.computebeta.func.V2.py

#step2.compute beta using reported:
python $rootdir/0_script/MRS.step2.computebeta.reported.py
```

```py
import pandas as pd
trait='BMI'
df1=pd.read_csv(trait+'.revise.lm005.csv',sep=',')
df1.shape
#(3513, 3)
df2=pd.read_csv(trait+'.revise.lmall.csv',sep=',')
df2.shape
#(3513, 3)
df3=pd.read_csv(trait+'.revise.lasso.csv',sep=',')
df3.shape
#(516411, 4)
dfx=pd.read_csv(trait+'.revise.reported.csv',sep=',')
dfx.shape
#3513,3

#id MRS kfold alpha
alphalist = df3['alpha'].drop_duplicates().tolist()
#147 alpha
for alphai in alphalist:
    print(alphai)
    df3i = df3[df3['alpha']==alphai]
    df3i = df3i[['id','kfold','MRS']]
    df3i.columns=['id','kfold','MRS_alpha_'+str(alphai)]
    print(df3i.shape)
    if alphai == alphalist[0]:
        df4 = df3i
    else:
        df4 = pd.merge(df4,df3i,on=['id','kfold'],how='inner')
    print(df4.shape)

df4.to_csv(trait+'.revise.lasso.allMRS.csv',sep=',',index=False)
#3513
df1 = df1[['id','kfold','MRS']]
df1.columns=['id','kfold','MRS_lm005']
df4 = pd.merge(df4,df1,on=['id','kfold'],how='inner')
#3513
df2 = df2[['id','kfold','MRS']]
df1.columns=['id','kfold','MRS_lmall']
df4 = pd.merge(df4,df1,on=['id','kfold'],how='inner')
#3513
dfx = dfx[['id','kfold','MRS']]
dfx.columns=['id','kfold','MRS_reported']
df4 = pd.merge(df4,dfx,on=['id','kfold'],how='inner')
#3513
phefile='/data/run/pheno.csv'
phe=pd.read_csv(phefile,sep=',')
phe['id']=phe.apply(lambda x: str(x['slide'])+'_'+x['array'],axis=1)
#3513
phe['id'].value_counts()
phe[phe['id']=='204087940132_R01C01']
phe.iloc[814,-1] = '204087940132_R01C01_1'
phe.iloc[1504,-1] = '204087940132_R01C01_2'
phe['id'].value_counts()
for phei in ['sbp','dbp','BMI']:
  meani=phe[phei].mean()
  phe[phei] = phe[phei].fillna(meani)

df4 = pd.merge(phe,df4,on=['id'],how='inner')
df4.to_csv(trait+'.revise.allMRS.csv',sep=',',index=False)
```
```sh
#MRS selection
#1st best lasso alpha
Rscript MRS.step3.analysis.func.R

#select best MRS in test
Rscript MRS.step4.compare.R
```
# compute MRS in STS1071
```sh
#1st extract cpgs from STS1071
Rscript /data1/yaning/work/PRSprofile/0_script/MRS.step5.extractfromSTS.R

#2st compute MRS in STS1071 and get file names=reg id  MRS target.resid    r2  pred
python /data1/yaning/work/PRSprofile/0_script/MRS.step6.computeMRSinSTS.py

#3st r2 in sts
Rscript /data1/yaning/work/PRSprofile/0_script/MRS.step6.r2inSTS.R
python /data1/yaning/work/PRSprofile/0_script/MRS.step6.r2inSTS.step2.R

```

```R
library('ggplot2')
library(tidyr)
library("tidyverse")
#library(ggthemes)
#library(RColorBrewer)
#library(colorspace)

setwd('E:\\04Zhangyaning\\PRS\\2_PRS_V3\\10_MRS\\MRS.revise')
mydata <- read.table('ALL.result.compare.V2.csv',sep=',',header=T)
mydata1 <- mydata[c('Trait', 'Method',  'TrainingR2','TrainingR2_ci005','TrainingR2_ci095')]
mydata1['population']='Tuning'
mydata2 <- mydata[c('Trait', 'Method',  'TestingR2','TestingR2_ci005','TestingR2_ci095')]
mydata2['population']='Testing'
mydata3 <- mydata[c('Trait', 'Method',  'ValidationR2','ValidationR2_ci005','ValidationR2_ci095')]
mydata3['population']='Validation'

names(mydata1) <- c('Trait', 'Method',  'R2','R2_CI5','R2_CI95','population')
names(mydata2) <- c('Trait', 'Method',  'R2','R2_CI5','R2_CI95','population')
names(mydata3) <- c('Trait', 'Method',  'R2','R2_CI5','R2_CI95','population')
newz <- rbind(mydata2,mydata3)
newz <- newz %>%  mutate(population = factor(population, levels = c("Tuning","Testing","Validation")))

newz <- newz %>%
  mutate(Method = factor(Method, levels = c("lasso", "lm_005", "lmall")))
levels(newz$Method)

#xx <- subset(newz,Trait=="BMI" & population=="Validation")

p <- ggplot(newz, aes(x = Method, y = R2, fill=Method)) + 
    geom_bar(color = 'black',position=position_dodge(), stat="identity") +
    scale_fill_manual(values = c('#C6E2FF','#8B8682','#8B0000'))+
    #geom_errorbar(aes(ymin=R2_CI5, ymax=R2_CI95), width=.4,    position=position_dodge(.9))+
    theme_bw() +
    theme( 
        #legend.position="none",
        axis.ticks = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.title.x=element_text(face="bold", size=35, color = "black", vjust = -0.5),  
        axis.title.y=element_text(face="bold", size=35, color = "black"),  
        axis.text.x = element_blank(), 
        axis.text.y=element_text(size=15, face="bold", color = "black")) +
    scale_y_continuous(limits = c(0, 0.2),breaks = seq(0,0.2,0.02))

p + facet_grid(population ~ Trait)+
    theme(
      strip.text.x = element_text(
        size = 20, color = "black", face = "bold"
        ), # 这里设置x轴方向的字体类型，
      strip.text.y = element_text(
        size = 20, color = "black", face = "bold"
        ) # 这里设置y轴方向的字体类型，
      ) +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave('E:\\04Zhangyaning\\PRS\\2_PRS_V3\\10_MRS\\MRS.revise\\testing.bestMRS.tiff', plot = last_plot(),
  scale = 1, width = 50, height = 25, units =c("cm"),
  dpi = 300, limitsize = TRUE)

setEPS()
postscript("E:\\04Zhangyaning\\PRS\\2_PRS_V3\\10_MRS\\MRS.revise\\testing.bestMRS.eps")
dev.off()
```

```sh
#3st MRS新的主图，mrsboxplot&OR&AUC for fat, HTN and OrH
python /data1/yaning/work/PRSprofile/0_script/MRS.step7.FigbestMRS.all.step1.py
$Rpwd /data1/yaning/work/PRSprofile/0_script/MRS.step7.FigbestMRS.all.step2.R
python /data1/yaning/work/PRSprofile/0_script/MRS.step7.FigbestMRS.all.AUC.V3.py

```
#计算r2&95%CI
```sh
Rpwd=/work/home/acd2j8na2s/soft/R/4.0.3/bin/Rscript
$Rpwd $rootdir/0_script//data1/yaning/work/PRSprofile/0_script/MRS.step8.RwithandwithoutsexandageforMRS.R
```
#95%CI AUC
```R
library(boot)
library(pROC)
# boot的使用方式很奇怪
get_auc <- function(data, ind, outcome, predictor){
  d = data[ind,] #这句必须加
  au <- as.numeric(auc(pROC::roc(d[,outcome], d[,predictor],quiet=T)))
  au
}
#HTN_definition
pop<-'STS1071'
for (target in c('OrH','fat','HTN_definition')) {
    print(target)
    phenofilefile <- paste('E:/04Zhangyaning/PRS/2_PRS_V3/10_MRS/MRS.revise/',pop,target,'.prob.csv',sep = "")
    dt = read.csv(phenofilefile, header = T, row.names = 1)
    get_auc(dt, outcome="target",predictor="M1_prob.1")
    set.seed(123)
    ba <- boot(dt, get_auc, R = 1000,
           outcome="target",predictor="M1_prob.1")
    print(ba)
    print(boot.ci(ba,conf = 0.95, type = "perc"))
}


```
#MRS for sex boxplot
```py
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np 
import os 
import seaborn as sns 
from scipy import stats
import statsmodels.api as sm
from  scipy.stats import chi2_contingency
import numpy as np
from scipy.stats import t
import statsmodels.formula.api as smf
import sys
from sklearn.preprocessing import StandardScaler
import matplotlib.font_manager as font_manager
import matplotlib as mpl
import random
import shutil
from functools import reduce
from statsmodels.formula.api import ols
from scipy import interp
import scipy.stats as st
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LinearRegression,LogisticRegression
from sklearn.model_selection import LeaveOneOut,StratifiedKFold,KFold
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.metrics import median_absolute_error, r2_score,confusion_matrix,explained_variance_score,mean_absolute_error
from statannotations.Annotator import Annotator
#pip install statannotations
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import KFold,StratifiedKFold

plt.switch_backend('agg')
plt.cla()
plt.close("all")
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300

def getPRSadjage(data,PRSname):
    X = data['age']
    y = data[PRSname]
    X = sm.add_constant(X)
    model = sm.OLS(y, X)
    results = model.fit()
    residuals = results.resid
    return(residuals)

font_size = 30
plt.rcParams['font.size'] = font_size
plt.rcParams['figure.figsize'] = (12,6)
mpl.rc('axes', lw=5)
workdir='/data1/yaning/work/PRSprofile/10_MRS/MRS.revise'
fig = plt.figure()
#subplot(numRows, numCols, plotNum)
axes = fig.subplots(1,2)

df=pd.read_csv('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/STS1071OrH.prob.csv')
dat1=df[['id','sex','MRS_BMI']]
dat2=df[['id','sex','MRS_HTN']]
dat1['MRS_BMIadj'] = getPRSadjage(df,'MRS_BMI')
dat2['MRS_HTNadj'] = getPRSadjage(df,'MRS_HTN')

ax1 = axes[0]
prsscorename='MRS_BMIadj'
c1, c0, c2 = sns.color_palette('Set1', 3)
dist1 = dat1[dat1['sex']==1][prsscorename]
dist0 = dat1[dat1['sex']==2][prsscorename]
mean1 = dist1.mean()
std1 = dist1.std()
mean0 = dist0.mean()
std0 = dist0.std()
print(mean0,std0)
print(mean1,std1)
stat_val01, p_val01 = stats.ttest_ind(dist0, dist1, equal_var=False)
p_val01 = "{:.2e}".format(p_val01)
print(p_val01)
my_pal = {1: "#f9766e", 2:"lightgreen"}
sns.boxplot(data=dat1, x='sex', y=prsscorename,ax=ax1,linewidth=5,palette=my_pal,width=0.8)
ax1.set_xticklabels(['Male','Female'],rotation = 0,fontsize = font_size)
pairs=[(1,2)]
annotator = Annotator(ax1, pairs, data=dat1, x='sex', y=prsscorename)
annotator.configure(test='Mann-Whitney', text_format='star',line_height=0.03,line_width=5)
annotator.apply_and_annotate()
ax1.set_ylabel('MRS_BMI',fontsize=font_size,fontweight='bold')
ax1.set_xlabel('Sex',fontsize=font_size,fontweight='bold')
ax1.set_ylim([-3,6])

ax2 = axes[1]
prsscorename='MRS_HTNadj'
c1, c0, c2 = sns.color_palette('Set1', 3)
dist1 = dat2[dat2['sex']==1][prsscorename]
dist0 = dat2[dat2['sex']==2][prsscorename]
mean1 = dist1.mean()
std1 = dist1.std()
mean0 = dist0.mean()
std0 = dist0.std()
print(mean0,std0)
print(mean1,std1)
stat_val01, p_val01 = stats.ttest_ind(dist0, dist1, equal_var=False)
p_val01 = "{:.2e}".format(p_val01)
print(p_val01)
my_pal = {1: "#f9766e", 2:"lightgreen"}
sns.boxplot(data=dat2, x='sex', y=prsscorename,ax=ax2,linewidth=5,palette=my_pal,width=0.8)
ax2.set_xticklabels(['Male','Female'],rotation = 0,fontsize = font_size)
pairs=[(1,2)]
annotator = Annotator(ax2, pairs, data=dat2, x='sex', y=prsscorename)
annotator.configure(test='Mann-Whitney', text_format='star',line_height=0.03,line_width=5)
annotator.apply_and_annotate()
ax2.set_ylabel('MRS_HTN',fontsize=font_size,fontweight='bold')
ax2.set_xlabel('Sex',fontsize=font_size,fontweight='bold')
ax2.set_ylim([-3,6])
plt.tight_layout()
plt.savefig('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/MRS.sexplot.png',bbox_inches='tight')
```

###  multi-omics R2+AUC
```sh
R $rootdir/0_script//data1/yaning/work/PRSprofile/0_script/Multi.omics.step1.Rwithandwithoutsexandage.R
python $rootdir/0_script//data1/yaning/work/PRSprofile/0_script/Multi.omics.step2.AUC.py
```
#95%CI AUC
```R
library(boot)
library(pROC)
# boot的使用方式很奇怪
get_auc <- function(data, ind, outcome, predictor){
  d = data[ind,] #这句必须加
  au <- as.numeric(auc(pROC::roc(d[,outcome], d[,predictor],quiet=T)))
  au
}
#HTN_definition
pop<-'STS1071'
for (target in c('OrH','fat','HTN_definition')) {
    print(target)
    phenofilefile <- paste('E:/04Zhangyaning/PRS/2_PRS_V3/11_MultiOmics/',pop,target,'.prob.csv',sep = "")
    dt = read.csv(phenofilefile, header = T, row.names = 1)
    get_auc(dt, outcome="target",predictor="M3_prob.1")
    set.seed(123)
    ba <- boot(dt, get_auc, R = 1000,
           outcome="target",predictor="M3_prob.1")
    print(ba)
    print(boot.ci(ba,conf = 0.95, type = "perc"))
}


```


```sh
#histogram
python /data1/yaning/work/PRSprofile/0_script/Multi.omics.step3.histogram.py HTN_definition HTN
python /data1/yaning/work/PRSprofile/0_script/Multi.omics.step3.histogram.py fat Obesity
python /data1/yaning/work/PRSprofile/0_script/Multi.omics.step3.histogram.py OrH OrH
```
```sh
cd /data1/yaning/work/PRSprofile/PGScatalog
###BMI
mkdir /data1/yaning/work/PRSprofile/PGScatalog/BMI.step1.download
cd /data1/yaning/work/PRSprofile/PGScatalog/BMI.step1.download
for pgsid in PGS000027 PGS000298 PGS000320 PGS000770 PGS000841 PGS000910 PGS000921 PGS002251 PGS002313 PGS002360 PGS002385 PGS002434 PGS002483 PGS002532 PGS002581 PGS002630 PGS002679 PGS002751 PGS002840 PGS002841 PGS002842 PGS002843 PGS002844 PGS002845 PGS002846 PGS002847 PGS002848 PGS002849 PGS002850 PGS002851 PGS002852 PGS002853 PGS002854 PGS002855 PGS002856 PGS002857 PGS002858 PGS002859 PGS003462 PGS003842 PGS003843 PGS003844 PGS003845 PGS003846 PGS003847 PGS003848 PGS003884 PGS003885 PGS003886 PGS003887 PGS003980 PGS003996 PGS004012 PGS004022 PGS004037 PGS004050 PGS004066 PGS004080 PGS004096 PGS004104 PGS004120 PGS004134 PGS004150 PGS004319 PGS004609 PGS004610 PGS004864 PGS004902 PGS004982 PGS004983 PGS004984 PGS004985 PGS004986 PGS004987 PGS004988 PGS004989 PGS004990 PGS004991 PGS004992 PGS004993 PGS004994
do
    wget https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${pgsid}/ScoringFiles/${pgsid}.txt.gz
done

for pgsid in PGS000027 PGS000298 PGS000320 PGS000770 PGS000841 PGS000910 PGS000921 PGS002251 PGS002313 PGS002360 PGS002385 PGS002434 PGS002483 PGS002532 PGS002581 PGS002630 PGS002679 PGS002751 PGS002840 PGS002841 PGS002842 PGS002843 PGS002844 PGS002845 PGS002846 PGS002847 PGS002848 PGS002849 PGS002850 PGS002851 PGS002852 PGS002853 PGS002854 PGS002855 PGS002856 PGS002857 PGS002858 PGS002859 PGS003462 PGS003842 PGS003843 PGS003844 PGS003845 PGS003846 PGS003847 PGS003848 PGS003884 PGS003885 PGS003886 PGS003887 PGS003980 PGS003996 PGS004012 PGS004022 PGS004037 PGS004050 PGS004066 PGS004080 PGS004096 PGS004104 PGS004120 PGS004134 PGS004150 PGS004319 PGS004609 PGS004610 PGS004864 PGS004902 PGS004982 PGS004983 PGS004984 PGS004985 PGS004986 PGS004987 PGS004988 PGS004989 PGS004990 PGS004991 PGS004992 PGS004993 PGS004994
do
gzip -d /data1/yaning/work/PRSprofile/PGScatalog/BMI.step1.download/${pgsid}.txt.gz
done

mkdir /data1/yaning/work/PRSprofile/PGScatalog/BMI.step2.remainSTS
python /data1/yaning/work/PRSprofile/0_script/PGScatalog.step2.remainSTS.py

#统计the number of 剩余items的各自remain SNP
ls /data1/yaning/work/PRSprofile/PGScatalog/BMI.step2.remainSTS > /data1/yaning/work/PRSprofile/PGScatalog/PGScatalog.BMI.step2.remainSTS.remainPGS.txt
python /data1/yaning/work/PRSprofile/0_script/PGScatalog.step2.remainSTS.remainSNP.py

mkdir /data1/yaning/work/PRSprofile/PGScatalog/BMI.step3.PRS
#####CASPMI
plinkpwd2='/data1/yaning/software/plink2/plink2'
for pgsid in `awk '{print $1}' /data1/yaning/work/PRSprofile/PGScatalog/PGScatalog.BMI.step2.remainSTS.remainPGS.txt | awk '{split($0,a,".");print a[1]}'`
do
    echo $pgsid
    $plinkpwd2 --bfile /data1/yaning/work/PRSprofile/1_Data/CASPMI --score /data1/yaning/work/PRSprofile/PGScatalog/BMI.step2.remainSTS/${pgsid}.txt 2 4 6 header list-variants cols=scoresums -out /data1/yaning/work/PRSprofile/PGScatalog/BMI.step3.PRS/${pgsid}
done
#####STS1071
mkdir /data1/yaning/work/PRSprofile/PGScatalog/BMI.step3.PRS/STS1071/
for pgsid in `awk '{print $1}' /data1/yaning/work/PRSprofile/PGScatalog/PGScatalog.BMI.step2.remainSTS.remainPGS.txt | awk '{split($0,a,".");print a[1]}'`
do
    echo $pgsid
    $plinkpwd2 --bfile /data1/yaning/work/PRSprofile/1_Data/STS1071 --score /data1/yaning/work/PRSprofile/PGScatalog/BMI.step2.remainSTS/${pgsid}.txt 2 4 6 header list-variants cols=scoresums -out /data1/yaning/work/PRSprofile/PGScatalog/BMI.step3.PRS/STS1071/${pgsid}
done
############################################################SBP
mkdir /data1/yaning/work/PRSprofile/PGScatalog/SBP.step1.download
cd /data1/yaning/work/PRSprofile/PGScatalog/SBP.step1.download
for pgsid in PGS000301 PGS000913 PGS002238 PGS002257 PGS002275 PGS002349 PGS002376 PGS002421 PGS002470 PGS002519 PGS002568 PGS002617 PGS002666 PGS002715 PGS002734 PGS002807 PGS003069 PGS003070 PGS003071 PGS003072 PGS003073 PGS003074 PGS003075 PGS003076 PGS003077 PGS003078 PGS003079 PGS003080 PGS003081 PGS003082 PGS003083 PGS003084 PGS003085 PGS003086 PGS003087 PGS003088 PGS003588 PGS003882 PGS003968 PGS003970 PGS003971 PGS004231 PGS004235 PGS004594 PGS004603 PGS005008 PGS005009 PGS005010 PGS005011 PGS005012 PGS005013 PGS005014 PGS005015 PGS005016 PGS005017 PGS005018 PGS005019 PGS005020
do
    wget https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${pgsid}/ScoringFiles/${pgsid}.txt.gz
done

for pgsid in PGS000301 PGS000913 PGS002238 PGS002257 PGS002275 PGS002349 PGS002376 PGS002421 PGS002470 PGS002519 PGS002568 PGS002617 PGS002666 PGS002715 PGS002734 PGS002807 PGS003069 PGS003070 PGS003071 PGS003072 PGS003073 PGS003074 PGS003075 PGS003076 PGS003077 PGS003078 PGS003079 PGS003080 PGS003081 PGS003082 PGS003083 PGS003084 PGS003085 PGS003086 PGS003087 PGS003088 PGS003588 PGS003882 PGS003968 PGS003970 PGS003971 PGS004231 PGS004235 PGS004594 PGS004603 PGS005008 PGS005009 PGS005010 PGS005011 PGS005012 PGS005013 PGS005014 PGS005015 PGS005016 PGS005017 PGS005018 PGS005019 PGS005020
do
gzip -d /data1/yaning/work/PRSprofile/PGScatalog/SBP.step1.download/${pgsid}.txt.gz
done

mkdir /data1/yaning/work/PRSprofile/PGScatalog/SBP.step2.remainSTS
python /data1/yaning/work/PRSprofile/0_script/PGScatalog.step2.remainSTS.py
#统计the number of 剩余items的各自remain SNP
ls /data1/yaning/work/PRSprofile/PGScatalog/SBP.step2.remainSTS > /data1/yaning/work/PRSprofile/PGScatalog/PGScatalog.SBP.step2.remainSTS.remainPGS.txt
python /data1/yaning/work/PRSprofile/0_script/PGScatalog.step2.remainSTS.remainSNP.py

mkdir /data1/yaning/work/PRSprofile/PGScatalog/SBP.step3.PRS
plinkpwd2='/data1/yaning/software/plink2/plink2'
for pgsid in `awk '{print $1}' /data1/yaning/work/PRSprofile/PGScatalog/PGScatalog.SBP.step2.remainSTS.remainPGS.txt | awk '{split($0,a,".");print a[1]}'`
do
    echo $pgsid
    $plinkpwd2 --bfile /data1/yaning/work/PRSprofile/1_Data/CASPMI --score /data1/yaning/work/PRSprofile/PGScatalog/SBP.step2.remainSTS/${pgsid}.txt 2 4 6 header list-variants cols=scoresums -out /data1/yaning/work/PRSprofile/PGScatalog/SBP.step3.PRS/${pgsid}
done
#####STS1071
mkdir /data1/yaning/work/PRSprofile/PGScatalog/SBP.step3.PRS/STS1071
plinkpwd2='/data1/yaning/software/plink2/plink2'
for pgsid in `awk '{print $1}' /data1/yaning/work/PRSprofile/PGScatalog/PGScatalog.SBP.step2.remainSTS.remainPGS.txt | awk '{split($0,a,".");print a[1]}'`
do
    echo $pgsid
    $plinkpwd2 --bfile /data1/yaning/work/PRSprofile/1_Data/STS1071 --score /data1/yaning/work/PRSprofile/PGScatalog/SBP.step2.remainSTS/${pgsid}.txt 2 4 6 header list-variants cols=scoresums -out /data1/yaning/work/PRSprofile/PGScatalog/SBP.step3.PRS/STS1071/${pgsid}
done
############################################################DBP
mkdir /data1/yaning/work/PRSprofile/PGScatalog/DBP.step1.download
cd /data1/yaning/work/PRSprofile/PGScatalog/DBP.step1.download
for pgsid in PGS000302 PGS000912 PGS002239 PGS002258 PGS002322 PGS002362 PGS002394 PGS002443 PGS002492 PGS002541 PGS002590 PGS002639 PGS002688 PGS002889 PGS002890 PGS002891 PGS002892 PGS002893 PGS002894 PGS002895 PGS002896 PGS002897 PGS002898 PGS002899 PGS002900 PGS002901 PGS002902 PGS002903 PGS002904 PGS002905 PGS002906 PGS002907 PGS002908 PGS003463 PGS003883 PGS003964 PGS004232 PGS004233 PGS004604
do
    wget https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${pgsid}/ScoringFiles/${pgsid}.txt.gz
done

for pgsid in PGS000302 PGS000912 PGS002239 PGS002258 PGS002322 PGS002362 PGS002394 PGS002443 PGS002492 PGS002541 PGS002590 PGS002639 PGS002688 PGS002889 PGS002890 PGS002891 PGS002892 PGS002893 PGS002894 PGS002895 PGS002896 PGS002897 PGS002898 PGS002899 PGS002900 PGS002901 PGS002902 PGS002903 PGS002904 PGS002905 PGS002906 PGS002907 PGS002908 PGS003463 PGS003883 PGS003964 PGS004232 PGS004233 PGS004604
do
gzip -d /data1/yaning/work/PRSprofile/PGScatalog/DBP.step1.download/${pgsid}.txt.gz
done

mkdir /data1/yaning/work/PRSprofile/PGScatalog/DBP.step2.remainSTS
python /data1/yaning/work/PRSprofile/0_script/PGScatalog.step2.remainSTS.py
#统计the number of 剩余items的各自remain SNP
ls /data1/yaning/work/PRSprofile/PGScatalog/DBP.step2.remainSTS > /data1/yaning/work/PRSprofile/PGScatalog/PGScatalog.DBP.step2.remainSTS.remainPGS.txt
python /data1/yaning/work/PRSprofile/0_script/PGScatalog.step2.remainSTS.remainSNP.py

mkdir /data1/yaning/work/PRSprofile/PGScatalog/DBP.step3.PRS
plinkpwd2='/data1/yaning/software/plink2/plink2'
for pgsid in `awk '{print $1}' /data1/yaning/work/PRSprofile/PGScatalog/PGScatalog.DBP.step2.remainSTS.remainPGS.txt | awk '{split($0,a,".");print a[1]}'`
do
    echo $pgsid
    $plinkpwd2 --bfile /data1/yaning/work/PRSprofile/1_Data/CASPMI --score /data1/yaning/work/PRSprofile/PGScatalog/DBP.step2.remainSTS/${pgsid}.txt 2 4 6 header list-variants cols=scoresums -out /data1/yaning/work/PRSprofile/PGScatalog/DBP.step3.PRS/${pgsid}
done
#####STS1071
mkdir /data1/yaning/work/PRSprofile/PGScatalog/DBP.step3.PRS/STS1071
plinkpwd2='/data1/yaning/software/plink2/plink2'
for pgsid in `awk '{print $1}' /data1/yaning/work/PRSprofile/PGScatalog/PGScatalog.DBP.step2.remainSTS.remainPGS.txt | awk '{split($0,a,".");print a[1]}'`
do
    echo $pgsid
    $plinkpwd2 --bfile /data1/yaning/work/PRSprofile/1_Data/STS1071 --score /data1/yaning/work/PRSprofile/PGScatalog/DBP.step2.remainSTS/${pgsid}.txt 2 4 6 header list-variants cols=scoresums -out /data1/yaning/work/PRSprofile/PGScatalog/DBP.step3.PRS/STS1071/${pgsid}
done
###############
mkdir /data1/yaning/work/PRSprofile/PGScatalog/step4.compare
#get file #"reg","id","PRS",'age','sex','target'
python /data1/yaning/work/PRSprofile/0_script/PGScatalog.step3.compare.py


for pgsid in `awk '{print $1}' /data1/yaning/work/PRSprofile/PGScatalog/PGScatalog.BMI.step2.remainSTS.remainPGS.txt | awk '{split($0,a,".");print a[1]}'`
do
    for pop in CASPMI STS1071 
    do
        Rscript /data1/yaning/work/PRSprofile/0_script/PGScatalog.step3.compare.R BMI $pgsid $pop
    done
done

for pgsid in `awk '{print $1}' /data1/yaning/work/PRSprofile/PGScatalog/PGScatalog.DBP.step2.remainSTS.remainPGS.txt | awk '{split($0,a,".");print a[1]}'`
do
    for pop in CASPMI STS1071 
    do
        Rscript /data1/yaning/work/PRSprofile/0_script/PGScatalog.step3.compare.R DBP $pgsid $pop
    done
done

for pgsid in `awk '{print $1}' /data1/yaning/work/PRSprofile/PGScatalog/PGScatalog.SBP.step2.remainSTS.remainPGS.txt | awk '{split($0,a,".");print a[1]}'`
do
    for pop in CASPMI STS1071 
    do
        Rscript /data1/yaning/work/PRSprofile/0_script/PGScatalog.step3.compare.R SBP $pgsid $pop
    done
done

for pop in CASPMI STS1071 
    do
        Rscript /data1/yaning/work/PRSprofile/0_script/PGScatalog.step3.compare.R BMI PRScsx $pop
        Rscript /data1/yaning/work/PRSprofile/0_script/PGScatalog.step3.compare.R DBP PRScsx $pop
        Rscript /data1/yaning/work/PRSprofile/0_script/PGScatalog.step3.compare.R SBP PRScsx $pop
    done

python /data1/yaning/work/PRSprofile/0_script/PGScatalog.step4.compare.py

```
我们观测到的男性更高的肥胖和高血压患病率更多的可能是因为表观遗传因素，可能的环境暴露探索
```sh
#健康生活方式的histplot
python /data1/yaning/work/PRSprofile/0_script/lifestyle.py

#代谢综合征的多个组分比如甘油三酯，空腹血糖，HDL，LDL等+lifestyle的6个组分的的性别差异
python /data1/yaning/work/PRSprofile/0_script/DNAm.MS.V3.py

#根据MRS预测值和实际值将样本分成前10%（低估），bottom（准确预测），后10%（高估），分析MRS高估的个体血糖显著升高，甘油三酯升高，高密度脂蛋白胆固醇降低。
python /data1/yaning/work/PRSprofile/0_script/DNAm.real.predict.py
```
在不同的cohort中男女患病率的比较
```sh
python /data1/yaning/work/PRSprofile/0_script/cohort.prevalence.sex.py
```

```py
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np 
import os 
import seaborn as sns 
from scipy import stats
from  scipy.stats import chi2_contingency
import numpy as np
from scipy.stats import t
import statsmodels.formula.api as smf
import sys
from sklearn.preprocessing import StandardScaler
import matplotlib.font_manager as font_manager
import matplotlib as mpl
import random
import shutil
from functools import reduce
from statsmodels.formula.api import ols
from scipy import interp
import scipy.stats as st
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LinearRegression,LogisticRegression
from sklearn.model_selection import LeaveOneOut,StratifiedKFold,KFold
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.metrics import median_absolute_error, r2_score,confusion_matrix,explained_variance_score,mean_absolute_error
from statannotations.Annotator import Annotator
#pip install statannotations
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import KFold,StratifiedKFold

plt.switch_backend('agg')
plt.cla()
plt.close("all")
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300

font_size = 30
plt.rcParams['font.size'] = font_size
plt.rcParams['figure.figsize'] = (8,6)
mpl.rc('axes', lw=5)
workdir='/data1/yaning/work/PRSprofile/10_MRS/MRS.revise'
fig, ax = plt.subplots()
dat=pd.read_csv('/data1/yaning/work/PRSprofile/8_bestPRS/STS1071OrH.prob.csv')
my_pal = {1: "orangered", 2: "lightgreen"}
dat['PRS for BMI'] = dat['PRS_BMI']
sns.boxplot(data=dat, x='sex', y='PRS for BMI',linewidth=5,palette=my_pal,width=0.8)
ax.set_xticklabels(['Male','Female'],rotation = 0,fontsize = font_size)
plt.tight_layout()
plt.savefig('/data1/yaning/work/PRSprofile/PRS_BMI.sexplot.png',bbox_inches='tight')

dist0 = dat[dat['sex']==1]['PRS for BMI']
dist1 = dat[dat['sex']==2]['PRS for BMI']
mean1 = dist1.mean()
std1 = dist1.std()
mean0 = dist0.mean()
std0 = dist0.std()
print(mean0,std0)
print(mean1,std1)
stat_val01, p_val01 = stats.ttest_ind(dist0, dist1, equal_var=False)
p_val01 = "{:.2e}".format(p_val01)
print(p_val01)

```
正态性检验
```R
phenofilefile <- paste('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/STS1071OrH.prob.csv',sep = "")
dt = read.csv(phenofilefile, header = T, row.names = 1)
ks.test(dt$PRS_BMI, "pnorm")
#D = 0.014736, p-value = 0.9742
#由于 p 值大于 0.05，我们接受原假设。我们有足够的证据表明样本数据来自正态分布
ks.test(dt$PRS_HTN, "pnorm")
#D = 0.02898, p-value = 0.3294
ks.test(dt$MRS_BMI, "pnorm")
#D = 0.018107, p-value = 0.8739
ks.test(dt$MRS_HTN, "pnorm")
#D = 0.022636, p-value = 0.6427

#[1] 2.716000e-13 1.042000e-05 1.333333e-02 4.200000e-01
p.adjust(
  c(1.27e-4,3.00e-42  ,0.01,0.42),  # P值列表,需要从小到大
  method ="bonferroni",n=4                       # FDR校正的方法
)
#[1] 5.08e-04 1.20e-41 4.00e-02 1.00e+00
```

cpg summary
```py
import pandas as pd 
import numpy as np
df1=pd.read_csv('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/selCPGsfrompbulished.txt',sep='\t')

trait='dbp'
trait2='DBP'
df2=pd.read_csv('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/'+trait+'.revise.lasso.best.summary.csv',sep=',')
df3=pd.read_csv('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/'+trait+'.revise.lm.summary005.csv',sep=',')
df4=pd.read_csv('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/'+trait+'.revise.lm.summaryall.csv',sep=',')
df11=df1[df1['Trait']==trait2]
df2list=df2['Features'].tolist()
df3list=df3['Unnamed: 0'].tolist()
df4list=df4['Unnamed: 0'].tolist()
df11['selectinlasso'] = df11.apply(lambda x: 'YES' if x['CpG'] in df2list else 'NO',axis=1)
df11['selectinlm005'] = df11.apply(lambda x: 'YES' if x['CpG'] in df3list else 'NO',axis=1)
df11['selectinlm(inNSPTandCAScohort)'] = df11.apply(lambda x: 'YES' if x['CpG'] in df4list else 'NO',axis=1)
df11.to_csv('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/'+trait+'.cpg.summary.csv',sep=',',index=False)

```

```R
library(bigsnpr)
library(bigreadr)
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(boot)
epicmapinfo="/data1/yaning/Data/STS/DNAm/epicmapinfo.RData"
#epicmapinfo: DNAm map position data file. Default: "/liufanGroup/public_data/gabDNAmData/epicmapinfo.RData" (EPIC array CpG map position hg19).
#DNAm CHR   BP_hg19
cpgfile='/data1/yaning/work/PRSprofile/allcps.txt'
cpg <- fread(cpgfile,sep='\t')
allcpg <- cpg$CpG

load(epicmapinfo)
ind.cpg  <- which(epicmapinfo$DNAm %in% allcpg)
epicmapinfo2 <- epicmapinfo[ind.cpg]
#313 cpgs= 234 for bmi + 44 for dbp + 74 for sbp
write.table (epicmapinfo2, file ='/data1/yaning/work/PRSprofile/epicmapinfo2', sep =",", row.names =FALSE, col.names =TRUE, quote =FALSE)
```
```py
import pandas as pd 
import numpy as np 
import os 
df1=pd.read_csv('/data1/yaning/work/PRSprofile/allcps.txt',sep='\t')
df2=pd.read_csv('/data1/yaning/work/PRSprofile/epicmapinfo2',sep=',')
df3 = pd.merge(df1,df2,how='left',left_on='CpG',right_on='DNAm')
df3.to_csv('/data1/yaning/work/PRSprofile/cpginsts.csv',sep=',',index=False)
```
```R
library(ChAMP)
library(data.table)
myimport <- champ.import(directory=system.file("extdata",package="ChAMPdata"))
myImport=myimport#包里的演示代码有个小细节错了，没有区分大小写，无伤大雅的
myfilter <- champ.filter(beta=myImport$beta,pd=myImport$pd,detP=myImport$detP,beadcount=myImport$beadcount)
str(hm450.manifest.hg19)

cpgfile='/data1/yaning/work/PRSprofile/file01.txt'
cpg <- fread(cpgfile,sep='\t')
ind.cpg  <- which(hm450.manifest.hg19$probeID %in% cpg$CpG)
x <- hm450.manifest.hg19[ind.cpg,c("probeID","CpG_chrm","CpG_beg","CpG_end",
          "gene_HGNC")]
#1579
write.table (x, file ='/data1/yaning/work/PRSprofile/allcpg.gene.txt', sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
```
```py
import pandas as pd 
import numpy as np 
df1=pd.read_csv('/data1/yaning/work/PRSprofile/file01.txt',sep='\t')
df2=pd.read_csv('/data1/yaning/work/PRSprofile/allcpg.gene.txt', sep ="\t")
df1[df1['CpG'].isin(df2['probeID'].tolist())==False].shape
#1579, cg03067296,cg00214628 have no gene anovargene;
df3 = pd.merge(df1,df2,left_on=['CpG'],right_on=['probeID'],how='left')
df3.to_csv('/data1/yaning/work/PRSprofile/cpginsts.gene.csv',sep=',',index=False)
```

#MRS for sex scatter plot
```py
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np 
import os 
import seaborn as sns 
from scipy import stats
import statsmodels.api as sm
from  scipy.stats import chi2_contingency
import numpy as np
from scipy.stats import t
import statsmodels.formula.api as smf
import sys
from sklearn.preprocessing import StandardScaler
import matplotlib.font_manager as font_manager
import matplotlib as mpl
import random
import shutil
from functools import reduce
from statsmodels.formula.api import ols
from scipy import interp
import scipy.stats as st
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LinearRegression,LogisticRegression
from sklearn.model_selection import LeaveOneOut,StratifiedKFold,KFold
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.metrics import median_absolute_error, r2_score,confusion_matrix,explained_variance_score,mean_absolute_error
from statannotations.Annotator import Annotator
#pip install statannotations
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import KFold,StratifiedKFold

plt.switch_backend('agg')
plt.cla()
plt.close("all")
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300
font_size = 30
plt.rcParams['font.size'] = font_size
plt.rcParams['figure.figsize'] = (8,6)
mpl.rc('axes', lw=5)
workdir='/data1/yaning/work/PRSprofile/10_MRS/MRS.revise'
fig, ax = plt.subplots()
dat=pd.read_csv('/data1/yaning/work/PRSprofile/8_bestPRS/STS1071OrH.prob.csv')
my_pal = {1: "orangered", 0: "lightgreen"}
dat['PRS for BMI'] = dat['PRS_BMI']
sns.boxplot(data=dat, x='sex', y='PRS for BMI',linewidth=5,palette=my_pal,width=0.8,hue='target',order=[1,2])
ax.set_xticklabels(['Male','Female'],rotation = 0,fontsize = font_size)
plt.tight_layout()
leg = plt.legend()
ax.get_legend().remove()
plt.savefig('/data1/yaning/work/PRSprofile/PRS_BMI.sexplot2.png',bbox_inches='tight')

fig, ax = plt.subplots()
dat=pd.read_csv('/data1/yaning/work/PRSprofile/8_bestPRS/STS1071OrH.prob.csv')
my_pal = {1: "orangered", 0: "lightgreen"}
dat['PRS for HTN'] = dat['PRS_HTN']
sns.boxplot(data=dat, x='sex', y='PRS for HTN',linewidth=5,palette=my_pal,width=0.8,hue='target',order=[1,2])
ax.set_xticklabels(['Male','Female'],rotation = 0,fontsize = font_size)
plt.tight_layout()
leg = plt.legend()
ax.get_legend().remove()
plt.savefig('/data1/yaning/work/PRSprofile/PRS_HTN.sexplot2.png',bbox_inches='tight')

fig, ax = plt.subplots()
dat=pd.read_csv('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/STS1071OrH.prob.csv')
my_pal = {1: "orangered", 0: "lightgreen"}
dat['MRS for BMI'] = dat['MRS_BMI']
sns.boxplot(data=dat, x='sex', y='MRS for BMI',linewidth=5,palette=my_pal,width=0.8,hue='target',order=[1,2])
ax.set_xticklabels(['Male','Female'],rotation = 0,fontsize = font_size)
plt.tight_layout()
leg = plt.legend()
ax.get_legend().remove()
plt.savefig('/data1/yaning/work/PRSprofile/MRS_BMI.sexplot2.png',bbox_inches='tight')
fig, ax = plt.subplots()
dat=pd.read_csv('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/STS1071OrH.prob.csv')
my_pal = {1: "orangered", 0: "lightgreen"}
dat['MRS for HTN'] = dat['MRS_HTN']
sns.boxplot(data=dat, x='sex', y='MRS for HTN',linewidth=5,palette=my_pal,width=0.8,hue='target',order=[1,2])
ax.set_xticklabels(['Male','Female'],rotation = 0,fontsize = font_size)
plt.tight_layout()
leg = plt.legend()
ax.get_legend().remove()
plt.savefig('/data1/yaning/work/PRSprofile/MRS_HTN.sexplot2.png',bbox_inches='tight')
```

#scatter plot
```py
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np 
import os 
import seaborn as sns 
from scipy import stats
import scipy
import statsmodels.api as sm
from  scipy.stats import chi2_contingency
import numpy as np
from scipy.stats import t
import statsmodels.formula.api as smf
import sys
from sklearn.preprocessing import StandardScaler
import matplotlib.font_manager as font_manager
import matplotlib as mpl
import random
import shutil
from functools import reduce
from statsmodels.formula.api import ols
from scipy import interp
import scipy.stats as st
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import LinearRegression,LogisticRegression
from sklearn.model_selection import LeaveOneOut,StratifiedKFold,KFold
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.metrics import median_absolute_error, r2_score,confusion_matrix,explained_variance_score,mean_absolute_error
from statannotations.Annotator import Annotator
#pip install statannotations
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import KFold,StratifiedKFold

plt.switch_backend('agg')
plt.cla()
plt.close("all")
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300
font_size = 30
plt.rcParams['font.size'] = font_size
plt.rcParams['figure.figsize'] = (8,6)
mpl.rc('axes', lw=5)
workdir='/data1/yaning/work/PRSprofile/10_MRS/MRS.revise'

df2 = pd.read_csv('/data1/yaning/work/PRSprofile/2_Phenotype/STS1071.PheforPaper.step2.csv',sep=',')
df2=df2[['id','BMI','dbp','sbp']]


fig, ax = plt.subplots()
dat=pd.read_csv('/data1/yaning/work/PRSprofile/8_bestPRS/STS1071OrH.prob.csv')
dat = pd.merge(dat,df2,on=['id'])
dat['PRS for BMI'] = dat['PRS_BMI']
#sns.set(style="whitegrid",font_scale=1.2)
my_pal = {1: "lightcoral", 2: "c"}
g=sns.lmplot(x='PRS_BMI', y='BMI', data=dat,
             hue='sex',palette=my_pal)
#g.fig.set_size_inches(10,8)
plt.text(-2, 35, 'y=1.00x+25.11 ',fontsize=15,color='red') 
plt.text(1, 17, 'y=0.93x+22.38',fontsize=15,color='royalblue') 
plt.savefig('/data1/yaning/work/PRSprofile/PRS_BMI.sexscatter.png',bbox_inches='tight')

fig, ax = plt.subplots()
dat=pd.read_csv('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/STS1071OrH.prob.csv')
dat = pd.merge(dat,df2,on=['id'])
dat['MRS for BMI'] = dat['MRS_BMI']
my_pal = {1: "lightcoral", 2: "c"}
g=sns.lmplot(x='MRS_BMI', y='BMI', data=dat,
             hue='sex',palette=my_pal)
plt.text(-2, 35, 'y=1.04x+24.99',fontsize=15,color='red') 
plt.text(1, 17, 'y=1.05x+22.56',fontsize=15,color='royalblue') 
plt.savefig('/data1/yaning/work/PRSprofile/MRS_BMI.sexscatter.png',bbox_inches='tight')

#################################################################################################
plt.switch_backend('agg')
plt.cla()
plt.close("all")
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300
font_size = 30
plt.rcParams['font.size'] = font_size
plt.rcParams['figure.figsize'] = (8,6)
mpl.rc('axes', lw=5)
workdir='/data1/yaning/work/PRSprofile/10_MRS/MRS.revise'

df2 = pd.read_csv('/data1/yaning/work/PRSprofile/2_Phenotype/STS1071.PheforPaper.step2.csv',sep=',')
df2=df2[['id','BMI','dbp','sbp']]
df2.columns=['id','BMI','DBP','SBP']

fig, ax = plt.subplots()
dat=pd.read_csv('/data1/yaning/work/PRSprofile/8_bestPRS/STS1071OrH.prob.csv')
dat = pd.merge(dat,df2,on=['id'])
my_pal = {1: "lightcoral", 2: "c"}
g=sns.lmplot(x='PRS_HTN', y='DBP', data=dat,
             hue='sex',palette=my_pal)

plt.text(-2, 110, 'y=1.81x+83.76',fontsize=15,color='red') 
plt.text(1, 75, 'y=1.95x+75.98',fontsize=15,color='royalblue') 
plt.savefig('/data1/yaning/work/PRSprofile/PRS_DBP.sexscatter.png',bbox_inches='tight')

fig, ax = plt.subplots()
dat=pd.read_csv('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/STS1071OrH.prob.csv')
dat = pd.merge(dat,df2,on=['id'])
my_pal = {1: "lightcoral", 2: "c"}
g=sns.lmplot(x='MRS_HTN', y='DBP', data=dat,
             hue='sex',palette=my_pal)
plt.text(-2, 110, 'y = 3.09x+82.63',fontsize=15,color='red') 
plt.text(1, 75, 'y = 2.50x+77.38',fontsize=15,color='royalblue') 
plt.savefig('/data1/yaning/work/PRSprofile/MRS_DBP.sexscatter.png',bbox_inches='tight')


#################################################################################################
plt.switch_backend('agg')
plt.cla()
plt.close("all")
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300
font_size = 30
plt.rcParams['font.size'] = font_size
plt.rcParams['figure.figsize'] = (8,6)
mpl.rc('axes', lw=5)
workdir='/data1/yaning/work/PRSprofile/10_MRS/MRS.revise'

df2 = pd.read_csv('/data1/yaning/work/PRSprofile/2_Phenotype/STS1071.PheforPaper.step2.csv',sep=',')
df2=df2[['id','BMI','dbp','sbp']]
df2.columns=['id','BMI','DBP','SBP']

fig, ax = plt.subplots()
#sns.set(style="whitegrid",font_scale=1.2)
dat=pd.read_csv('/data1/yaning/work/PRSprofile/8_bestPRS/STS1071OrH.prob.csv')
dat = pd.merge(dat,df2,on=['id'])
my_pal = {1: "lightcoral", 2: "c"}
g=sns.lmplot(x='PRS_HTN', y='SBP', data=dat,
             hue='sex',palette=my_pal)
plt.text(-2, 160, 'y=2.47x+125.38',fontsize=15,color='red') 
plt.text(1, 100, 'y=1.97x+112.55',fontsize=15,color='royalblue') 
plt.savefig('/data1/yaning/work/PRSprofile/PRS_SBP.sexscatter.png',bbox_inches='tight')

fig, ax = plt.subplots()
dat=pd.read_csv('/data1/yaning/work/PRSprofile/10_MRS/MRS.revise/STS1071OrH.prob.csv')
dat = pd.merge(dat,df2,on=['id'])
my_pal = {1: "lightcoral", 2: "c"}
g=sns.lmplot(x='MRS_HTN', y='SBP', data=dat,
             hue='sex',palette=my_pal)
#add regression equation to plot
plt.text(-2, 160, 'y = 4.34x+123.80',fontsize=15,color='red') 
plt.text(1, 100, 'y = 3.85x+114.61',fontsize=15,color='royalblue') 
plt.savefig('/data1/yaning/work/PRSprofile/MRS_SBP.sexscatter.png',bbox_inches='tight')
```
```py
##########################################################################################
df2 = pd.read_csv('/data1/yaning/work/PRSprofile/2_Phenotype/STS1071.PheforPaper.step2.csv',sep=',')
df2=df2[['id','BMI','dbp','sbp']]
df2.columns=['id','BMI','DBP','SBP']

fig, ax = plt.subplots()
dat=pd.read_csv('/data1/yaning/work/PRSprofile/8_bestPRS/STS1071OrH.prob.csv')
dat = pd.merge(dat,df2,on=['id'])
my_pal = {1: "orangered", 2: "lightgreen"}
#g=sns.lmplot(x='MRS_HTN', y='SBP', data=dat, hue='sex')
#create regplot
df1=dat[dat['sex']==1]
g1 = sns.regplot(data=df1, x=df1.PRS_HTN , y=df1.SBP,color='orangered')
#calculate slope and intercept of regression equation
slope, intercept, r, p, sterr = scipy.stats.linregress(x= g1.get_lines()[0].get_xdata(),
                                                       y=g1.get_lines()[0].get_ydata())
#add regression equation to plot
plt.text(-2, 160, 'y = ' + str(round(intercept,2)) + '+' + str(round(slope,2)) + 'x',fontsize=15,color='red')
y=2.47x+125.38 
#create regplot
fig, ax = plt.subplots()
dat=pd.read_csv('/data1/yaning/work/PRSprofile/8_bestPRS/STS1071OrH.prob.csv')
dat = pd.merge(dat,df2,on=['id'])
my_pal = {1: "orangered", 2: "lightgreen"}
df2=dat[dat['sex']==2]
g2 = sns.regplot(data=df2, x=df2.PRS_HTN , y=df2.SBP,color='blue')
#calculate slope and intercept of regression equation
slope2, intercept2, r2, p2, sterr2 = scipy.stats.linregress(x= g2.get_lines()[0].get_xdata(),
                                                       y=g2.get_lines()[0].get_ydata())
y=1.97x+112.55
```


