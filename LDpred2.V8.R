library(bigsnpr)
library(bigreadr)
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(boot)


args=commandArgs(T)
rootdir <- args[1]
db <- args[2]
trait <- args[3]
pop <- args[4]

#rootdir <- '/work/home/acd2j8na2s/Work/PRSprofile'
#db <- 'UKBB'
#trait <- 'BMI'
#pop <- 'EUR'

restdir <- paste(rootdir,'7_single_ancestry_PRS',trait,'ldpred2',db,sep = "/")
lddir <- paste(rootdir,'7_single_ancestry_PRS',trait,'LDmatrices',pop,sep = "/")
genfile <- paste(rootdir,'5_allGWASandData',trait,'plink',sep = "/")
plinkfile <- paste(genfile,'ALL',sep = "/")
pcafile <- paste(genfile,'ALL.eigenvec',sep = "/")
phenofilefile <- paste(rootdir,'7_single_ancestry_PRS/ALL.phe.txt',sep = "/")
gwasfile <- paste(rootdir,'5_allGWASandData',trait,db,paste(trait,'gwas.QC.CAS', sep='.'),sep = "/")
ldchr <- paste(lddir,paste(pop,'LD_chr', sep=''),sep = "/")
w_hm3 <- paste(rootdir,'0_script','w_hm3.snplist',sep = "/")
#'/work/home/acd2j8na2s/Work/PRSprofile/0_script/w_hm3.snplist'


if(trait=='SBP') {
    affection <- 'sbp'
}
if(trait=='DBP') {
    affection <- 'dbp'
}
if(trait=='BMI') {
    affection <- 'BMI'
}

#rootdir <- '/work/home/acrt8eogli/Work/PRSprofile'
#restdir <- '/work/home/acrt8eogli/Work/PRSprofile/7_single_ancestry_PRS/BMI/ldpred2/BBJ'
#lddir <- '/work/home/acrt8eogli/Work/PRSprofile/7_single_ancestry_PRS/BMI/LDmatrices/EAS'
#plinkfile <- '/work/home/acrt8eogli/Work/PRSprofile/5_allGWASandData/BMI/plink/ALL'
#gwasfile <- '/work/home/acrt8eogli/Work/PRSprofile/5_allGWASandData/BMI/BBJ/BMI.gwas.QC.CAS'
#phenofilefile <- '/work/home/acrt8eogli/Work/PRSprofile/7_single_ancestry_PRS/ALL.phe.txt'
#pcafile <- '/work/home/acrt8eogli/Work/PRSprofile/5_allGWASandData/BMI/plink/ALL.eigenvec'
#affection <- 'BMI'
#ldchr <- '/work/home/acrt8eogli/Work/PRSprofile/7_single_ancestry_PRS/BMI/LDmatrices/EAS/EASLD_chr'
############################################### rds data #########################################################################
#snp_readBed(paste(plinkfile,'bed',sep='.'))
#datatest <- snp_attach(paste(plinkfile,'rds',sep='.'))
#datatest <- snp_attach('/work/home/acrt8eogli/Work/PRSprofile/7_single_ancestry_PRS/BMI/ldpred2/BBJ/maptest.rds')
datatest <- snp_attach(paste(rootdir,'7_single_ancestry_PRS',trait,'ldpred2',db,'maptest.rds',sep = "/"))

Gtest <- datatest$genotypes
Gtest <- snp_fastImputeSimple(Gtest, method = "mean2")

map_test <- dplyr::transmute(datatest$map,chr = as.integer(chromosome), pos = physical.pos,
                        a0 = allele2, a1 = allele1)  # reversed somehow..
CHRtest <- as.integer(map_test$chr)

fam.order <- as.data.table(datatest$fam)
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
  file = "/work/home/acd2j8na2s/Work/PRSprofile/1_Data/STS2030.fam", 
  sep = " ", header = F)
sts1071fam <- read.table(
  file = "/work/home/acd2j8na2s/Work/PRSprofile/1_Data/STS1071.fam", 
  sep = " ", header = F)
caspmifam <- read.table(
  file = "/work/home/acd2j8na2s/Work/PRSprofile/1_Data/CASPMI.fam", 
  sep = " ", header = F)

#ix <- which(datatest$fam$family.ID %in% c(stsfam$V1))
#ind.val <- ix
#ind.test <- setdiff(rows_along(Gtest), ix)
ind.val  <- which(datatest$fam$family.ID %in% c(sts2030fam$V1))
ind.test  <- which(datatest$fam$family.ID %in% c(caspmifam$V1))
ind.extend  <- which(datatest$fam$family.ID %in% c(sts1071fam$V1))

################################################# Information for the variants provided in the LD reference #########################################################################################
## Information for the variants provided in the LD reference
#map_ldref <- readRDS(paste(lddir,"map.rds",sep='/'))
#dataref <- snp_attach('/work/home/acrt8eogli/Work/PRSprofile/7_single_ancestry_PRS/BMI/ldpred2/BBJ/mapref.rds')
map_ldref <- readRDS(paste(rootdir,'7_single_ancestry_PRS',trait,'ldpred2',db,'mapref.rds',sep = "/"))
dim(map_ldref)
######################################################### summary statistics #####################################################
# Read external summary statistics
sumstats <- fread2(gwasfile, na.strings = "NULL",
                   select = c("SNP","CHR", "BP", "OA", "EA", "Beta", "se",'N'),
                   col.names = c("rsid","chr", "pos", "a0", "a1", "beta", "beta_se","n_eff"))

# Read w_hm3
w_hm3df <- read.csv(w_hm3,sep='\t',header=T)
#1217311 SNPs

# Filter out tets and ref and w_hm3 SNPs
in_test <- vctrs::vec_in(sumstats[, c("chr", "pos")], map_test[, c("chr", "pos")])
sumstats <- sumstats[in_test, ]
dim(sumstats)
#2707939
in_test2 <- vctrs::vec_in(sumstats[, c("rsid")], w_hm3df[, c("SNP")])
sumstats <- sumstats[in_test2, ]
dim(sumstats)
#607362

info_snp <- snp_match(sumstats, map_ldref)
tmp <- tempfile(tmpdir = "/work/home/acd2j8na2s/Work/PRSprofile/tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# calculate LD
CHR <- as.integer(map_ldref$chr)
for (chr in 1:22) {
    cat(chr, ".. ", sep = "")
    # Extract SNPs that are included in the chromosome
    ## indices in 'info_snp'
    #info_snpex = info_snp[info_snp['rsid']=='rs1418707',]
    ind.chr <- which(info_snp$chr == chr)
    #rs1418707 in 231187,
    ## indices in 'Gref'
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    #rs1418707 in 1053006
    ## G indices in 'corr'
    ind.chr3 <- match(ind.chr2, which(CHR == chr))
    print('ok1')
    # read corr
    corr0 <- readRDS(paste0(ldchr, chr, ".rds"))
    #207560 SNPs
    corr0 <- corr0[ind.chr3, ind.chr3]
    print('ok2')
    if (chr == 1) {
      df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
        ld <- Matrix::colSums(corr0^2)
        print('ok3')
        # compact = TRUE: Stores only x, but all (even the zero ones) from first to last being not 0
        # Stores all (i, x) for x != 0
        corr <- as_SFBM(corr0,tmp, compact = TRUE)
        print('ok4')
    } else {
      df_beta <- rbind(df_beta, info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
        ld <- c(ld, Matrix::colSums(corr0^2))
        print('ok3')
        corr$add_columns(corr0, nrow(corr))
        print('ok4')
    }
}

# Heritability estimation of LD score regression
(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(df_beta),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = 15,blocks = NULL)))

h2_est <- ldsc[["h2"]]
########################################################################################################################
# LDpred2-inf
  #beta_inf <- snp_ldpred2_inf(corr, df_beta, h2_est)

# LDpred2-grid
  #(h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4))
  #h2_seq <- 0.1681 for BMI in BBJ
  (h2_seq <- round(h2_est * c(1), 4))
  #(p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
  #0.001,0.003,0.01,0.03,0.1,0.3,1
  (p_seq <- c(0.003,0.01,0.03,0.1,0.3,1))
  #Disease odds ratio by an increase of one standard deviation of GPS
  #RFs: clinical risk factors
  (params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))

  beta_grid <- snp_ldpred2_grid(corr, df_beta, params)
  params$sparsity <- colMeans(beta_grid == 0)

  #bigparallelr::set_blas_ncores(NCORES)
  pred_grid <- big_prodMat(Gtest, beta_grid, ind.row = ind.val,
                           ind.col = df_beta[["_NUM_ID_"]])
  #params$score <- big_univLogReg(as_FBM(pred_grid), y[ind.val])$score

  params$score <- apply(pred_grid, 2, function(x) {
  if (all(is.na(x))) return(NA)
  #summary(lm(y[ind.val] ~ x))$coef["x", 3]
  summary(lm(y[ind.val] ~ x))$r.squared
  # summary(glm(y[ind.val] ~ x, family = "binomial"))$coef["x", 3]
})

rsq <- function(data, indices){
  d <-  data[indices,]
  formula <- "target~PRS"
  fit <- lm(formula, data = d)
  return(summary(fit)$r.square)
 }

params$score_ci005 <- apply(pred_grid, 2, function(x) {
  if (all(is.na(x))) return(NA)
  reg.dat <- data.frame(
    target = y[ind.val],
    PRS = x,stringsAsFactors = FALSE)
  #reg.dat <- y[ind.val]
  #reg.dat$PRS <- x
  reg.datresults <- boot(data=reg.dat, statistic=rsq, R=1000)
  reg.datci=boot.ci(reg.datresults,type="norm")
  ci1r2 = c(reg.datci$normal)[2]
  })

params$score_ci095 <- apply(pred_grid, 2, function(x) {
  if (all(is.na(x))) return(NA)
  reg.dat <- data.frame(
    target = y[ind.val],
    PRS = x,stringsAsFactors = FALSE)
  #reg.dat <- y[ind.val]
  #reg.dat$PRS <- x
  reg.datresults <- boot(data=reg.dat, statistic=rsq, R=1000)
  reg.datci=boot.ci(reg.datresults,type="norm")
  ci2r2 = c(reg.datci$normal)[3]
  })

write.table (params, file =paste(restdir,'ldpred2grid.csv',sep='/'), sep =",", row.names =FALSE, col.names =TRUE, quote =FALSE)

  best_beta_grid_nosp <- params %>%
    mutate(id = row_number()) %>%
    filter(!sparse) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_grid[, .]

  best_beta_grid_sp <- params %>%
    mutate(id = row_number()) %>%
    filter(sparse) %>%
    arrange(desc(score)) %>%
    slice(1) %>%
    pull(id) %>%
    beta_grid[, .]

  # compute predictions for test set
  betas <- cbind(best_beta_grid_nosp, best_beta_grid_sp)
  pred_test <- big_prodMat(Gtest, betas, ind.row = ind.test,
                           ind.col = df_beta[["_NUM_ID_"]])
  pred_val <- big_prodMat(Gtest, betas, ind.row = ind.val,
                           ind.col = df_beta[["_NUM_ID_"]])
  pred_extend <- big_prodMat(Gtest, betas, ind.row = ind.extend,
                           ind.col = df_beta[["_NUM_ID_"]])



p <- ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  facet_wrap(~ sparse, labeller = label_both) +
  labs(y = "r.squared", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"),
    axis.ticks = element_blank(), 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(),
    axis.title.x=element_text(face="bold", size=35, color = "black", vjust = -0.5),  
    axis.title.y=element_text(face="bold", size=35, color = "black"),  
    axis.text.x=element_text(size=15, face="bold", color = "black", vjust=0.1), 
    axis.text.y=element_text(size=15, face="bold", color = "black"))

pdf(paste(restdir,'ldpred2grid.pdf',sep='/'),width = 10,height = 5)
p
dev.off()

# save results
#SNP list以及对应的EA,NonEA也需要保存
ind.col = df_beta[["_NUM_ID_"]]
savenewgwasinfo <- map_test[ind.col,]
savenewgwasinfo$best_beta_grid_nosp <- best_beta_grid_nosp
savenewgwasinfo$best_beta_grid_sp <- best_beta_grid_sp
saveRDS(savenewgwasinfo, paste(restdir,'ldpred2.savenewgwasinfo.rds',sep='/'))

#id,pred需要保存
dfval = data.frame(id = datatest$fam[ind.val,'sample.ID'], ldpred2_beta_grid_nosp = pred_val[,1], 
  ldpred2_beta_grid_sp = pred_val[,2])
dftest = data.frame(id = datatest$fam[ind.test,'sample.ID'], ldpred2_beta_grid_nosp = pred_test[,1], 
  ldpred2_beta_grid_sp = pred_test[,2])
dfextend = data.frame(id = datatest$fam[ind.extend,'sample.ID'], ldpred2_beta_grid_nosp = pred_extend[,1], 
  ldpred2_beta_grid_sp = pred_extend[,2])

z <- rbind(dfval, dftest,dfextend)
write.table (z, file =paste(restdir,'ldpred2.ALL.pred.csv',sep='/'), sep =",", row.names =FALSE, col.names =TRUE, quote =FALSE)
