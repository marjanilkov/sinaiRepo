# # Check the PC1 and PC2 against all phenotypes/covariates to see if there are 
# # any confounding variables that affect the outcome of the clustering which 
# # shows two very robust clusters 
# 
# rm(list = ls())
# library(WGCNA)
# library(Biobase)
# setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/hippoBlood_AgeSub/anew2/PC/")
# 
# # load the data to be checked
# blood = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/hippoBlood_AgeSub/old/20221017_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_log2CPM.Whole_blood.sex_rin_adj.RDS")
# blood_pheno = pData(blood)
# blood_expr = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/hippoBlood_AgeSub/anew/20221027_data/Whole_blood.sex_rin_ischTime_adj_std.RDS")
# 
# blood_expr = as.data.frame(t(blood_expr))
# # 
# # blood_pheno = merge(ME, blood_pheno, by.x = "row.names", by.y = "row.names")
# # rownames(blood_pheno) = blood_pheno$Row.names
# # blood_comb = merge(blood_pheno, blood_expr, by.x = "row.names", by.y = "row.names")
# 
# #plot the PCA
# blood.pca = prcomp(blood_expr,
#                    center = F,
#                    scale. = F)
# saveRDS(blood.pca, "blood.pca.ischTime.RDS")
# 
# blood.pca = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/hippoBlood_AgeSub/anew/20221031_PCvsCovariates/before_correction_for_SMCENTER/blood.pca.RDS")
# # cluster_colors = vector( "integer" , 56200 )
# # blood.pca.wgcna = moduleEigengenes(blood_expr, cluster_colors,
# #                                    nPC = 2, 
# #                                    softPower = 1,
# #                                    scale = F,
# #                                    verbose = 10, indent = 0)
# # saveRDS(blood.pca.wgcna, "blood.pca.wgcna.RDS")


rm(list = ls())
library(ggfortify)
library(cluster)
library(Biobase)
setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/20221114_GTEx_ADsub_biomark/20221128_PCA/")

# load the data to be checked
datExpr = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/20221114_GTEx_ADsub_biomark/20221124_data/Whole_Blood_adj_std.RDS")
clinical = read.delim("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/hippoBlood_AgeSub/anew2/20221116_data/GTEx_Analysis_v8_Annotations_clinical_covariates.txt")
datExpr = datExpr[,1:(ncol(datExpr)-8)]
# we need the cases where there are no NAs
clinical = clinical[complete.cases(clinical),]
# also SMCENTER has 10 rows with "", not NA but just an empty holder. We need to
# take care of that as well
clinical = clinical[!clinical$SMCENTER=="",]

# we need to add character to the DTHHRDY variable since it has values 0-4 which
# the linear model lm function later on takes it as numerical while in reality 
# this is a categorization variable which states the Hardy death scale
clinical$DTHHRDY = paste("H",clinical$DTHHRDY, sep = "")

# We will now extract the data we need
phenotype = "Whole_Blood"

df = clinical[clinical$SMTSD==phenotype, ]
#datExpr = gtex[,df$SAMPID]
pheno = clinical[clinical$SAMPID %in% df$SAMPID,]

#plot the PCA
#datExpr = as.data.frame(t(datExpr))

datExpr.pca = prcomp(datExpr,
                   center = F,
                   scale. = F)
#saveRDS(datExpr.pca, "Brain_HippocampusNOADJ.pca.RDS")



library(ggfortify)
library(cluster)


# fix some variables 
# blood_pheno$SEX = as.character(blood_pheno$SEX)
# blood_pheno$DTHHRDY = ifelse(is.na(blood_pheno$DTHHRDY), 5, blood_pheno$DTHHRDY)
# blood_pheno$DTHHRDY = as.character(blood_pheno$DTHHRDY)
tmp = c("SMCENTER", "SEX", "AGE", "DTHHRDY", "SMRIN", "SMTSISCH", "SMRRNART", "SMNTERRT")
#tmp_num = c("SMRIN", "SMTSISCH", "SMEXNCRT", "SMRRNART", "SMNTERRT")

tmp1 = "DTHHRDY"

#for (tmp1 in tmp){
  print(tmp1)
  autoplot(datExpr.pca, 
         data = pheno,
         colour = tmp1,
         size = 3)+
  ggsave(
    paste(tmp1,".png",sep = ""),
    plot = last_plot(),
    scale = 1,
    width = 250,
    height = 200,
    units = c( "mm"),
    dpi = 300,
    limitsize = TRUE)
#}

  blood_pheno$SMCENTER = ifelse(blood_pheno$SMCENTER == "B1", 1, blood_pheno$SMCENTER)
  blood_pheno$SMCENTER = ifelse(blood_pheno$SMCENTER == "C1", 2, blood_pheno$SMCENTER)
  blood_pheno$SMCENTER = ifelse(blood_pheno$SMCENTER == "D1", 3, blood_pheno$SMCENTER)
  cor(y$PC1,as.integer(blood_pheno$SMCENTER),method = "pearson")
  plot( blood_pheno$SMCENTER,y$PC1)



autoplot(blood.pca, 
         data = blood_pheno,
         colour = 'SMCENTER',
         size = 4)

lmSMTSISCH = lm(ME0 ~ SMTSISCH, data = blood_comb)
summary(lmSMTSISCH)
plot(blood_comb$ME0, blood_comb$SMTSISCH )

lmSMRIN = lm(ME0 ~ SMRIN, data = blood_comb)
summary(lmSMCENTER)
plot(blood_comb$ME0, blood_comb$SMRIN )

lmSMEXNCRT = lm(ME0 ~ SMEXNCRT, data = blood_comb)
summary(lmSMEXNCRT)
plot(blood_comb$ME0, blood_comb$SMEXNCRT )

lmSMRRNART = lm(ME0 ~ SMRRNART, data = blood_comb)
summary(lmSMRRNART)
plot(blood_comb$ME0, blood_comb$SMRRNART )

lmSMNTERRT = lm(ME0 ~ SMNTERRT, data = blood_comb)
summary(lmSMNTERRT)
plot(blood_comb$ME0, blood_comb$SMNTERRT )

lmAGE = lm(ME0 ~ AGE, data = blood_comb)
summary(lmAGE)
plot(blood_comb$ME0, blood_comb$AGE )


