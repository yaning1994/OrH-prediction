import pandas as pd 
import numpy as np 


for trait in ['BMI','DBP','SBP']:
    for db in ['BBJ','UKBB']:
        gwasfile='/work/home/acd2j8na2s/Work/PRSprofile/5_allGWASandData/'+trait+'/'+db+'/'+trait+'.gwas.QC.CAS'
        gwasdf = pd.read_csv(gwasfile,sep='\t')
        gwasdf = gwasdf[['SNP', 'CHR', 'EA', 'OA', 'Beta', 'se', 'N']]
        gwasdf.columns=['rsid','chr','a1','a0','beta','beta_se','n_eff']
        gwasdf.to_csv('/work/home/acd2j8na2s/Work/PRSprofile/7_single_ancestry_PRS/'+trait+'/PROSPER/'+db+'.GWAS.txt',sep='\t',index=False)

    