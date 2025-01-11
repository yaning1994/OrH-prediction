library(bigsnpr)
library(bigreadr)
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(boot)
library(tidyverse)

NCORES <- 30  # TO MODIFY
#args=commandArgs(T)
#rootdir <- args[1]
#db <- args[2]
#trait <- args[3]
#pop <- args[4]
rootdir <- '/work/home/acd2j8na2s/Work/PRSprofile'
restdir <- '/work/home/acd2j8na2s/Work/PRSprofile/test'
lddir <- '/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/BMI/LDmatrices/EAS'
plinkfile <- '/work/home/acd2j8na2s/Work/PRSprofile/test/testsample'
gwasfile <- '/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/BMI/BBJ/BMI.gwas.QC.CAS'
phenofilefile <- '/work/home/acd2j8na2s/Work/PRSprofile/test/test.phe.txt'
affection <- 'BMI'
refplinkfile <- '/work/home/acd2j8na2s/1000G/EAS_1kg'
pcafile <- '/work/home/acd2j8na2s/Work/PRSprofile/test/testsample.eigenvec'
#########################################################################
#snp_readBed(paste(refplinkfile,'bed',sep='.'))
dataref <- snp_attach(paste(refplinkfile,'rds',sep='.'))
map_ref <- dplyr::transmute(dataref$map,chr = as.integer(chromosome), pos = physical.pos,
                        a0 = allele2, a1 = allele1)  # reversed somehow..
#snp_readBed(paste(plinkfile,'bed',sep='.'))
datatest <- snp_attach(paste(plinkfile,'rds',sep='.'))
map_test <- dplyr::transmute(datatest$map,chr = as.integer(chromosome), pos = physical.pos,
                        a0 = allele2, a1 = allele1)  # reversed somehow..
in_ref <- vctrs::vec_in(map_ref[, c("chr", "pos")], map_test[, c("chr", "pos")])
X1 <- map_ref[, c("chr", "pos","a0",'a1')]
X1['beta']<-1
X2 <- map_test
X3 <- snp_match(X1, X2)
#Indices of the columns (SNPs) to keep
if (file.exists("/work/home/acd2j8na2s/Work/PRSprofile/test/mapref.bk"))
{
    file.remove("/work/home/acd2j8na2s/Work/PRSprofile/test/mapref.bk")
}
if (file.exists("/work/home/acd2j8na2s/Work/PRSprofile/test/mapref.rds"))
{
    file.remove("/work/home/acd2j8na2s/Work/PRSprofile/test/mapref.rds")
}

if (file.exists("/work/home/acd2j8na2s/Work/PRSprofile/test/maptest.bk"))
{
    file.remove("/work/home/acd2j8na2s/Work/PRSprofile/test/maptest.bk")
}
if (file.exists("/work/home/acd2j8na2s/Work/PRSprofile/test/maptest.rds"))
{
    file.remove("/work/home/acd2j8na2s/Work/PRSprofile/test/maptest.rds")
}

newdatarefile <- snp_subset(dataref,ind.col = X3$`_NUM_ID_.ss`,backingfile='/work/home/acd2j8na2s/Work/PRSprofile/test/mapref')
newdatatestile <- snp_subset(datatest,ind.col = X3$`_NUM_ID_`,backingfile='/work/home/acd2j8na2s/Work/PRSprofile/test/maptest')

newdataref <- snp_attach(newdatarefile)
newdatatest <- snp_attach(newdatatestile)
###################################################################################
Gref <- newdataref$genotypes
Gref <- snp_fastImputeSimple(Gref, method = "mean2", ncores = nb_cores())
map_ref <- dplyr::transmute(newdataref$map,chr = as.integer(chromosome), pos = physical.pos,
                        a0 = allele2, a1 = allele1)  # reversed somehow..
CHRref <- as.integer(newdataref$map$chromosome)
POSref <- newdataref$map$physical.pos


Gtest <- newdatatest$genotypes
Gtest <- snp_fastImputeSimple(Gtest, method = "mean2", ncores = nb_cores())
map_test <- dplyr::transmute(newdatatest$map,chr = as.integer(chromosome), pos = physical.pos,
                        a0 = allele2, a1 = allele1)  # reversed somehow..
CHRtest <- as.integer(newdatatest$map$chromosome)
POStest <- newdatatest$map$physical.pos

fam.order <- as.data.table(newdatatest$fam)
# We assume the fam order is the same across different chromosomes
# Rename fam order
setnames(fam.order,
        c("family.ID", "sample.ID"),
        c("FID", "IID"))

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
y <- pheno[fam.order, on = c("FID", "IID")]$target
#2030 STS2030 sampe index in merged as ind.train
sts2030fam <- read.table(
  file = "/work/home/acd2j8na2s/Work/PRSprofile/test/STS2030keepsample.txt", 
  sep = "\t", header = F)

sts1071fam <- read.table(
  file = "/work/home/acd2j8na2s/Work/PRSprofile/test/STS1071keepsample.txt", 
  sep = "\t", header = F)

caspmifam <- read.table(
  file = "/work/home/acd2j8na2s/Work/PRSprofile/test/CASPMIkeepsample.txt", 
  sep = "\t", header = F)

#正式的sep需要修改为空格
ind.val  <- which(datatest$fam$family.ID %in% c(sts2030fam$V1))
ind.test  <- which(datatest$fam$family.ID %in% c(caspmifam$V1))
ind.extend  <- which(datatest$fam$family.ID %in% c(sts1071fam$V1))
#ind.val <- ix
#ind.test <- setdiff(rows_along(Gtest), ix)
################################################# Information for the variants provided in the LD reference #########################################################################################
######################################################### summary statistics #####################################################
# Read external summary statistics
sumstats <- fread2(gwasfile, na.strings = "NULL",
                   select = c("SNP","CHR", "BP", "OA", "EA", "Beta", "se",'N','P'),
                   col.names = c("rsid","chr", "pos", "a0", "a1", "beta", "beta_se","n_eff",'p'))

in_test <- vctrs::vec_in(sumstats[, c("chr", "pos")], map_test[, c("chr", "pos")])
sumstats <- sumstats[in_test, ]
in_test <- vctrs::vec_in(sumstats[, c("chr", "pos")], map_ref[, c("chr", "pos")])
sumstats <- sumstats[in_test, ]

# Here, you also want to restrict to the variants present
# in your test data as well. For this, you can use something like
info_snpref <- snp_match(sumstats, map_ref)
# Rename column "ind.mapref" to "_NUM_ID_"
info_snpref <- rename(info_snpref, c("ind.mapref" = "_NUM_ID_"))

info_snptest <- snp_match(sumstats, map_test)
# Rename column "ind.mapref" to "_NUM_ID_"
info_snptest <- rename(info_snptest, c("ind.maptest" = "_NUM_ID_"))
############################################################################
### Clumping in reference genofile
# beta and lpval need to have the same length as ncol(G), CHR and POS
# -> one solution is to use missing values and use the 'exclude' parameter
lpvalref <- rep(NA, ncol(Gref))
lpvalref[info_snpref$`ind.mapref`] <- -log10(info_snpref$p)
all_keep <- snp_grid_clumping(Gref, CHRref, POSref, lpS = lpvalref,  
  grid.thr.r2 = c(0.01),
  grid.base.size = c(50), exclude = which(is.na(lpvalref)), ncores = NCORES)
attr(all_keep, "grid")
#SNP indices that are kept 
#for each chromosome, for each set of variants resulting from clumping and for each p-value threshold,
#这个index是Gref中的index，与Gtest一定一致
### Thresholding and make PRS in testgenofile
betatest <- rep(0, ncol(Gtest))
betatest[info_snptest$`ind.maptest`] <- info_snptest$beta

lpvaltest <- rep(NA, ncol(Gtest))
lpvaltest[info_snptest$`ind.maptest`] <- -log10(info_snptest$p)

if (file.exists(paste(restdir,'scores.rds.bk',sep='/'))) 
{
  #Delete file if it exists
  file.remove(paste(restdir,'scores.rds.rds',sep='/'))
  file.remove(paste(restdir,'scores.rds.bk',sep='/'))
}

#snp_grid_PRS-An FBM (matrix on disk) that stores the C+T scores for all parameters of the grid (and for each chromosome separately)
multi_PRS <- snp_grid_PRS(
      Gtest, all_keep, betas = betatest, lpS = lpvaltest, ind.row=ind.val,
      grid.lpS.thr = seq_log(0.1, 0.9999 * max(lpvaltest, na.rm = TRUE), 1),
      backingfile = paste(restdir,'scores.rds',sep='/'),
      ncores = NCORES
    )

str(multi_PRS)

### Stacking
final_mod <- snp_grid_stacking(multi_PRS, y[ind.val],  K = 5, ncores = NCORES)
mod <- final_mod$mod
summary(mod)
new_beta <- final_mod$beta.G

length(ind <- which(new_beta != 0))  # 491,240
summary(new_beta)
summary(new_beta[which(sign(new_beta * betatest) < 0)])

pred_sct_val <- final_mod$intercept + big_prodVec(Gtest, new_beta[ind], ind.row = ind.val, ind.col = ind)
pred_sct_test <- final_mod$intercept + big_prodVec(Gtest, new_beta[ind], ind.row = ind.test, ind.col = ind)
pred_sct_extend <- final_mod$intercept + big_prodVec(Gtest, new_beta[ind], ind.row = ind.extend, ind.col = ind)

### Best C+T
grid2 <- attr(all_keep, "grid") %>%
    mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), id = row_number()) %>%
    unnest(cols = "thr.lp")
s <- nrow(grid2)
grid2$R2 <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.val) {
    # Sum over all chromosomes, for the same C+T parameters
    single_PRS <- rowSums(X[, ind + s * (0:21)]) ## replace by 0:21 in real data
    reg.dat <- data.frame("target" = y.val, "PRS" = single_PRS)
    reg.formula <- "target ~ PRS"
    summary(lm(reg.formula, dat=reg.dat))$r.squared
  }, ind = 1:s, s = s, y.val = y[ind.val],
  a.combine = 'c', block.size = 1, ncores = NCORES)

return_two_value <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.val) {
    # Sum over all chromosomes, for the same C+T parameters
    library(boot)
    rsq <- function(data, indices){
        d <-  data[indices,]
        formula <- "target~PRS"
        fit <- lm(formula, data = d)
        return(summary(fit)$r.square)
    }
    single_PRS <- rowSums(X[, ind + s * (0:21)]) ## replace by 0:21 in real data
    reg.dat <- data.frame("target" = y.val, "PRS" = single_PRS)
    reg.formula <- "target ~ PRS"
    reg.datresults <- boot(data=reg.dat, statistic=rsq, R=1000)
    reg.datci=boot.ci(reg.datresults,type="norm")
    ci2r1 = c(reg.datci$normal)[2]
    ci2r2 = c(reg.datci$normal)[3]
    out <- list(one=ci2r1, two=ci2r2)
    return(out)
  }, ind = 1:s, s = s, y.val = y[ind.val],
  a.combine = 'c', block.size = 1, ncores = NCORES)
grid2$R2_ci005 <- return_two_value$one
grid2$R2_ci095 <- return_two_value$two
#并行工作每个都在一个干净的 R session 中运行，因此您必须加载 boot，rsq在每个 worker 中
##############################################################################################
max_prs <- grid2 %>% arrange(desc(R2)) %>% slice(1:10) %>% print() %>% slice(1)
ind.keep <- unlist(map(all_keep, max_prs$id))
pred_clumping_val <- snp_PRS(Gtest, betatest[ind.keep], ind.test = ind.val, ind.keep = ind.keep,
                           lpS.keep = lpvaltest[ind.keep], thr.list = max_prs$thr.lp)
pred_clumping_test <- snp_PRS(Gtest, betatest[ind.keep], ind.test = ind.test, ind.keep = ind.keep,
                           lpS.keep = lpvaltest[ind.keep], thr.list = max_prs$thr.lp)
pred_clumping_extend <- snp_PRS(Gtest, betatest[ind.keep], ind.test = ind.extend, ind.keep = ind.keep,
                           lpS.keep = lpvaltest[ind.keep], thr.list = max_prs$thr.lp)
# Save results
#grid
write.table (grid2, file =paste(restdir,'sctgrid.csv',sep='/'), sep =",", row.names =FALSE, col.names =TRUE, quote =FALSE)
#beta
#id,pred需要保存
dfval = data.frame(id = datatest$fam[ind.val,'sample.ID'], ct = pred_clumping_val, sct = pred_sct_val)
dftest = data.frame(id = datatest$fam[ind.test,'sample.ID'], ct = pred_clumping_test, sct = pred_sct_test)
dfextend = data.frame(id = datatest$fam[ind.extend,'sample.ID'], ct = pred_clumping_extend, sct = pred_sct_extend)
z <- rbind(dfval, dftest,dfextend)
write.table (z, file =paste(restdir,'ALL.pred.csv',sep='/'), sep =",", row.names =FALSE, col.names =TRUE, quote =FALSE)

#SNP list以及对应的EA,NonEA也需要保存
#clumpinginfo and sctinfo
savenewgwasinfo <- map_test
savenewgwasinfo$beta <- betatest
savenewgwasinfo$best_betact <- betatest
ix2 <- setdiff(rows_along(map_test), ind.keep)
savenewgwasinfo$best_betact[ix2] <- NA
savenewgwasinfo$best_betasct <- new_beta
saveRDS(savenewgwasinfo, paste(restdir,'savenewgwasinfo.rds',sep='/'))

