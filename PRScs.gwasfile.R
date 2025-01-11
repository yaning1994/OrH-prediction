library(dplyr)
args = commandArgs(trailingOnly=TRUE)


#a1 = '/thinker/storage/org/liufanGroup/zhangyn/FeasibilityAssessment/8_baseDataProcess' 
#a2 = 'UKBB.BMI.maf001.info08.gwas.alighned'

a1 = args[1]
a2 = args[2]


z.in <- paste0(a1, "/", a2)
z <- read.table(z.in, header = T,sep='\t')
z2 <- z[, c("SNP", "EA", "OA", "Beta", "P")]
colnames(z2) <- c("SNP", "A1", "A2", "BETA",  "P")
write.table(z2, paste0(a1, "/", a2,".PRScs.txt"), col.names = T, row.names = F, sep = "\t", quote = F)

#SNP          A1   A2   BETA      P
