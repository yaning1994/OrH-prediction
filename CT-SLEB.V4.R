##############################Data and Software Preparation
#install the CTSLEB package
#install.packages("devtools")
library(devtools)
## Loading required package: usethis
#install_github("andrewhaoyu/CTSLEB")
library(CTSLEB) #No
library(data.table)
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:data.table':
## 
##     between, first, last
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
library(caret)
## Loading required package: lattice
## Loading required package: ggplot2
library(SuperLearner)
## Loading required package: nnls
## Super Learner
## Version: 2.0-26
## Package created on 2019-10-27
#install.packages("ranger")
library(ranger) #No
library(glmnet)
## Loading required package: Matrix
## Loaded glmnet 4.1
library(boot)
#################################################################

args=commandArgs(T)
trait <- args[1]
#trait <- 'BMI'
if(trait=='SBP') {
    affection <- 'sbp'
}
if(trait=='DBP') {
    affection <- 'dbp'
}
if(trait=='BMI') {
    affection <- 'BMI'
}
#Specify the directory to the data folder
data_dir = "/work/home/acd2j8na2s/Work/PRSprofile"
#Specify the directory for the summary statistics
#EUR_sumstats_file <- '/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/BMI/UKBB/BMI.gwas.QC.CAS' 
#  reference population
EUR_sumstats_file <- paste(data_dir,'5_allGWASandData',trait,'UKBB',paste(trait,'.gwas.QC.CAS',sep=''),sep='/')

#AFR_sumstats_file <- '/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/BMI/BBJ/BMI.gwas.QC.CAS' 
#  target population
AFR_sumstats_file <- paste(data_dir,'5_allGWASandData',trait,'BBJ',paste(trait,'.gwas.QC.CAS',sep=''),sep='/')

#Specify the directory to the reference data for clumping purpose
EUR_ref_plinkfile <- paste(data_dir,'5_allGWASandData','allsnp.EUR',sep='/')
AFR_ref_plinkfile <- paste(data_dir,'5_allGWASandData','allsnp.EAS',sep='/')
#Specify the genotype data with subjects for tuning and validation
AFR_test_plinkfile <- paste(data_dir,'5_allGWASandData',trait,'plink','ALL',sep='/')
#Specify the directory to PLINK1.9 and PLINK2.0
plink19_exec <- "/work/home/acd2j8na2s/software/plink"
plink2_exec <- "/work/home/acd2j8na2s/software/plink2"
#Specify the directory to the result directory, must have '/' at last
out_dir = paste(data_dir,'7_single_ancestry_PRS',trait,'CTSLEB/',sep = "/")
outprefix <- "CTSLEB"
#Phenotypes 
phenofilefile <- paste(data_dir,'7_single_ancestry_PRS','ALL.phe.txt',sep='/')
pcafile <- paste(data_dir,'5_allGWASandData',trait,'plink','ALL.eigenvec',sep='/')
#################################
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
fam.order <- fread(paste(AFR_test_plinkfile,".fam",sep=''))
fam.order <- fam.order %>% select("V1", "V2") 
# We assume the fam order is the same across different chromosomes
# Rename fam order
setnames(fam.order,
        c("V1", "V2"),
        c("FID", "IID"))
y <- pheno[fam.order, on = c("FID", "IID")]$target
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
y_tune <- y[ind.val]
#Phenotypes for validation the target PRS 
y_vad <- y[ind.test]
#Phenotypes for test the target PRS 
y_test <- y[ind.extend]

##############################Step 1: Two-dimensional Clumping and Thresholding (CT)
#get SNP set with corresponding tuning parameters for estimating the covariance matrix for the prior distribution

#load data from the EUR and the target population
sum_EUR <- fread(EUR_sumstats_file,header=T)
sum_AFR <- fread(AFR_sumstats_file,header=T)
sum_EUR$SNP2 <- sum_EUR$SNP
sum_EUR <- sum_EUR %>% select("CHR", "SNP", 'BP','EA','Beta','se','P','SNP2') 
names(sum_EUR) <- c("CHR", "SNP", 'BP','A1',"BETA",'SE','P','rs_id')
head(sum_EUR)
##    CHR                      SNP       BP A1      BETA          SE      P
## 1:  22  rs55926024:16054740:A:G 16054740  G  0.003215 0.004537115 0.4786
##          rs_id
## 1:   rs2844885
sum_AFR$SNP2 <- sum_AFR$SNP
sum_AFR <- sum_AFR %>% select("CHR", "SNP", 'BP','EA','Beta','se','P','SNP2') 
names(sum_AFR) <- c("CHR", "SNP", 'BP','A1',"BETA",'SE','P','rs_id')
head(sum_AFR)
##    CHR                      SNP       BP A1      BETA         SE      P
## 1:  22  rs55926024:16054740:A:G 16054740  G -0.013930 0.01157938 0.2290
##         rs_id
## 1:  rs2844885
#Prepare the parameters for PLINK clumping purpose
#wc_base_vec = c(50, 100, 200, 500),
#r2_vec = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
#pthres = c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0),


PRS_farm <- SetParamsFarm(plink19_exec = plink19_exec,
                          plink2_exec = plink2_exec,
                          wc_base_vec = c(50, 100, 200, 500),
                          r2_vec = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95),
                          pthres = c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0),
                          threads = 2,
                          mem = 8000
     )


#'plink_list' not found untill now
#compute all prs by dimCT for all samples(tuning,validation,and tetsing)
#return 5 GlobalEnv: prs_mat,plink_list,snp_list,write_list,sum_com
prs_mat <- dimCT(results_dir = out_dir,
                 sum_target = sum_AFR,
                 sum_ref = sum_EUR,
                 ref_plink = EUR_ref_plinkfile,
                 target_plink = AFR_ref_plinkfile,
                 test_target_plink = AFR_test_plinkfile,
                 out_prefix = outprefix,
                 params_farm = PRS_farm)

str(plink_list)
#plink_list$scores, p_values, unique_infor, q_range, score_file, p_value_file, q_range_file
str(snp_list)
str(write_list)
str(sum_com)

dim(prs_mat)
#4092 10
names(prs_mat)
# [1] "#FID"                                          
# [2] "IID"                                           
# [3] "clump_r2_0.8_ws_62.5_p_other_5e-08_p_tar_5e-08"
# [4] "clump_r2_0.5_ws_100_p_other_5e-08_p_tar_5e-08" 
# [5] "clump_r2_0.8_ws_62.5_p_other_5e-08_p_tar_5e-07"
# [6] "clump_r2_0.5_ws_100_p_other_5e-08_p_tar_5e-07" 
# [7] "clump_r2_0.8_ws_62.5_p_other_5e-07_p_tar_5e-08"
# [8] "clump_r2_0.5_ws_100_p_other_5e-07_p_tar_5e-08" 
# [9] "clump_r2_0.8_ws_62.5_p_other_5e-07_p_tar_5e-07"
#[10] "clump_r2_0.5_ws_100_p_other_5e-07_p_tar_5e-07"


#prs_mat contains all samples. We will next use ind.val of these samples and calculate a set of tuning parameters
#Phenotypes for tuning the target PRS is located in “data/y_tuning.txt”. 
#The first two columns of prs_tune are the family id and individual ids. The target PRSs starts from the third column
prs_tune <- prs_mat[ind.val,]
n.total.prs <- length(pthres)^2*length(r2_vec)*length(wc_base_vec)
prs_r2_vec_test <- rep(0,n.total.prs)
prs_r2_CI5_vec_test <- rep(0,n.total.prs)
prs_r2_CI95_vec_test <- rep(0,n.total.prs)
rsq <- function(data, indices){
  d <-  data[indices,]
  formula <- "target~PRS"
  fit <- lm(formula, data = d)
  return(summary(fit)$r.square)
 }

for(p_ind in 1:n.total.prs){
  model <- lm(y_tune~prs_tune[,(2+p_ind)])
  rsqdata = data.frame(target = y_tune, PRS = prs_tune[,(2+p_ind)])
  prs_r2_vec_test[p_ind] <- summary(model)$r.square
  reg.datresults <- boot(data=rsqdata, statistic=rsq, R=1000)
  reg.datci=boot.ci(reg.datresults,type="norm")
  prs_r2_CI5_vec_test[p_ind] <- c(reg.datci$normal)[2]
  prs_r2_CI95_vec_test[p_ind] <- c(reg.datci$normal)[3]
}

max_ind <- which.max(prs_r2_vec_test)
print(colnames(prs_tune)[max_ind+2])
## [1] "clump_r2_0.01_ws_10000_p_other_5e-08_p_tar_0.05"
#it’s using clumping r2-cutoff at 0.01, window size at 10000kb, SNPs with p_EUR < 5E-08 or p_target < 0.05. 

#save grid in ind.val
gridcol1 <- colnames(prs_tune)[-c(1:2)]
grid = data.frame(PRSname = gridcol1, score = prs_r2_vec_test, score_ci005=prs_r2_CI5_vec_test,score_ci095=prs_r2_CI95_vec_test)
write.table (grid, file =paste(out_dir,'CTSLEB.grid.csv',sep=''), sep =",", row.names =FALSE, col.names =TRUE, quote =FALSE)

##############################Step 2: Empirical-Bayes (EB) Estimation of Effect Sizes
#find the best snp set, Actually, it's confirmed p in EAS GWAS,p in EUR GWAS,r2,ws
best_snps <- colnames(prs_tune)[max_ind+2]
#[1] "clump_r2_0.5_ws_100_p_other_5e-08_p_tar_5e-08"

#calculate eb effect using EB coefficients
#plink_list: List of plink data.frame and files from PreparePlinkFiles()
#params_farm: List of plink parameters produced from SetParamsFarm()
prs_mat_eb <- CalculateEBEffectSize(bfile = AFR_test_plinkfile ,
                                    snp_ind = best_snps,
                                    plink_list = plink_list, 
                                    out_prefix = outprefix,
                                    results_dir = out_dir,
                                    params_farm = PRS_farm)

dim(prs_mat_eb)
#[1] 4092 18
names(prs_mat_eb)
write.table (prs_mat_eb, file =paste(out_dir,'CTBLEB.prs_mat_eb.csv',sep=''), sep =",", row.names =FALSE, col.names =TRUE, quote =FALSE)
# [1] "#FID"                                                    
# [2] "IID"                                                     
# [3] "clump_r2_0.8_ws_62.5_EB_target_p_other_5e-08_p_tar_5e-08"
# [4] "clump_r2_0.8_ws_62.5_EB_ref_p_other_5e-08_p_tar_5e-08"   
# [5] "clump_r2_0.5_ws_100_EB_target_p_other_5e-08_p_tar_5e-08" 
# [6] "clump_r2_0.5_ws_100_EB_ref_p_other_5e-08_p_tar_5e-08"    
# [7] "clump_r2_0.8_ws_62.5_EB_target_p_other_5e-08_p_tar_5e-07"
# [8] "clump_r2_0.8_ws_62.5_EB_ref_p_other_5e-08_p_tar_5e-07"   
# [9] "clump_r2_0.5_ws_100_EB_target_p_other_5e-08_p_tar_5e-07" 
#[10] "clump_r2_0.5_ws_100_EB_ref_p_other_5e-08_p_tar_5e-07"    
#[11] "clump_r2_0.8_ws_62.5_EB_target_p_other_5e-07_p_tar_5e-08"
#[12] "clump_r2_0.8_ws_62.5_EB_ref_p_other_5e-07_p_tar_5e-08"   
#[13] "clump_r2_0.5_ws_100_EB_target_p_other_5e-07_p_tar_5e-08" 
#[14] "clump_r2_0.5_ws_100_EB_ref_p_other_5e-07_p_tar_5e-08"    
#[15] "clump_r2_0.8_ws_62.5_EB_target_p_other_5e-07_p_tar_5e-07"
#[16] "clump_r2_0.8_ws_62.5_EB_ref_p_other_5e-07_p_tar_5e-07"   
#[17] "clump_r2_0.5_ws_100_EB_target_p_other_5e-07_p_tar_5e-07" 
#[18] "clump_r2_0.5_ws_100_EB_ref_p_other_5e-07_p_tar_5e-07"
##############################Step 3: Super Learning
#Get the SNP set for estimating the covariance matrix of prior distribution in empirical Bayes method
prs_tune_sl <- prs_mat_eb[which(prs_mat_eb$IID %in% c(sts2030fam$V1)),c(-1,-2)]
prs_valid_sl <- prs_mat_eb[which(prs_mat_eb$IID %in% c(caspmifam$V1)),c(-1,-2)]
prs_test_sl <- prs_mat_eb[which(prs_mat_eb$IID %in% c(sts1071fam$V1)),c(-1,-2)]
dim(prs_tune_sl)

#Next train the super-learning model
SL.libray <- c(
  "SL.glmnet",
  "SL.ridge",
  "SL.nnet"
  #"SL.bayesglm"
  #"SL.stepAIC"
  #"SL.xgboost"
  #"SL.randomForest",
  #"SL.ksvm",
  #"SL.bartMachine",
  #"SL.kernelKnn",
  #"SL.rpartPrune",
  #"SL.lm"
  #"SL.mean"
)

sl <- SuperLearner(Y = y_tune, 
                   X = prs_tune_sl, 
                   family = gaussian(),
                   SL.library = SL.libray)

#Predict the outcome using the independent test dataset and evaluate the CT-SLEB PRS performance.
filnalr2 <- function(y_pred,y,rsq){
  model <- lm(y~y_pred)
  rsqdata = data.frame(target = y, PRS = y_pred)
  r2_ctsleb <- summary(model)$r.square
  reg.datresults <- boot(data=rsqdata, statistic=rsq, R=1000)
  reg.datci=boot.ci(reg.datresults,type="norm")
  r2_ctsleb_CI5 <- c(reg.datci$normal)[2]
  r2_ctsleb_CI95 <- c(reg.datci$normal)[3]
  c(r2_ctsleb,r2_ctsleb_CI5,r2_ctsleb_CI95)
}
y_tune_pred <- c(predict(sl, prs_tune_sl, onlySL = TRUE)[[1]])
y_valid_pred <- c(predict(sl, prs_valid_sl, onlySL = TRUE)[[1]])
y_test_pred <- c(predict(sl, prs_test_sl, onlySL = TRUE)[[1]])

tune_filnalr2 <- filnalr2(y_tune_pred,y_tune,rsq)
val_filnalr2 <- filnalr2(y_valid_pred,y_vad,rsq)
test_filnalr2 <- filnalr2(y_test_pred,y_test,rsq)

#id,pred需要保存
dfval = data.frame(id = c(sts2030fam$V1), CTSLEB = y_tune_pred)
dftest = data.frame(id = c(caspmifam$V1), CTSLEB = y_valid_pred)
dfextend = data.frame(id = c(sts1071fam$V1), CTSLEB = y_test_pred)
z <- rbind(dfval, dftest,dfextend)
write.table (z, file =paste(out_dir,'CTBLEB.ALL.pred.csv',sep=''), sep =",", row.names =FALSE, col.names =TRUE, quote =FALSE)

#r2需要保存
dfvalr = data.frame(R2 = tune_filnalr2[1],R2_ci005 = tune_filnalr2[2], R2_ci095 = tune_filnalr2[3])
dftestr = data.frame(R2 = val_filnalr2[1],R2_ci005 = val_filnalr2[2], R2_ci095 = val_filnalr2[3])
dfextendr = data.frame(R2 = test_filnalr2[1],R2_ci005 = test_filnalr2[2], R2_ci095 = test_filnalr2[3])
zr <- rbind(dfvalr, dftestr,dfextendr)
write.table (zr, file =paste(out_dir,'CTBLEB.filanR2.csv',sep=''), sep =",", row.names =FALSE, col.names =TRUE, quote =FALSE)
