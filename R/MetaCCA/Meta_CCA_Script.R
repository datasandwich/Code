library(metaCCA)

# Edit directory path
setwd('/Volumes/Google Drive/My Drive/PROJECT/Data/Aging Biomarkers/European/Jupyter_Exports/Individual Analysis')

#Edit the file path
S_XY_full = read.csv('S_XY_full/S_XY_full0.csv',row.names=1)
S_XY_filtered = read.csv('S_XY_filtered/S_XY_filtered0.csv',row.names=1)

levels(S_XY_full[,1]) = levels(S_XY_full_study1[,1])
levels(S_XY_full[,2]) = levels(S_XY_full_study1[,2])

levels(S_XY_filtered[,1]) = levels(S_XY_full_study1[,1])
levels(S_XY_filtered[,2]) = levels(S_XY_full_study1[,2])

S_YY_E = estimateSyy(S_XY_full)

N1 = lengths(S_XY_full['allele_0'], use.names=FALSE)

metaCCA_res1 = metaCcaGp( nr_studies = 1,
							S_XY = list(S_XY_filtered),
							std_info = c(0),
							S_YY = list(S_YY_E),
							N = c(N1) )
							
print(head(metaCCA_res1[1:2]), digits = 2)

# Edit the file name
write.csv(metaCCA_res1[1:2],'Results/Results_Chrom_1.csv')