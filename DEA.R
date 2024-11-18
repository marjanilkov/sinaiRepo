rm(list = ls())

options(stringsAsFactors=F)
library(limma)
################################################################################
# Variables set by the user
beta=1
dataMod = "meth"
p.val_cutoff =0.05
diffExpr = "DEM"
################################################################################
wkdir = paste("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew6/initial_cluster_out_test/",dataMod,".beta",beta,"/", sep = "");
setwd(wkdir)

# Read in the data
if (dataMod=="prot"){modSize = 9209; FC_cutoff = 1.2; datExpr <- readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/msbb.AMPAD.AD.proteomics.CDR_adjusted.uniq_genes.BM36.covcorr.ryan.std.RDS")}
if (dataMod=="mRNA"){modSize = 23201; FC_cutoff = 1.5; datExpr <- readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/BM_36.RNA.CDR_age_sex_adjusted.std.RDS")}
if (dataMod=="meth"){modSize = 270; FC_cutoff = 1.1; datExpr <- readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/BM_36.MSBB_DMR270_mean_exp.CDR_age_sex_adjusted.std.RDS")}

datExpr = datExpr[complete.cases(datExpr),]
classification = readRDS(paste("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew6/initial_cluster_out_test/",dataMod,".beta",beta,"/",dataMod,".subtypes.WINA.beta",beta,".RDS", sep = ""))

clusters = unique(classification$cluster)
datExpr = datExpr[,classification$SynapseId]

logCPM = datExpr
design = model.matrix(~0+cluster,data=classification) 

fit = lmFit(datExpr,design) #fit all groups together!
contr <- makeContrasts(C1_VS_ctrl = clusterC1 - clustercontrol, 
                       A_VS_ctrl = clusterA - clustercontrol, 
                       levels = colnames(coef(fit)))

vfit <- contrasts.fit(fit, contr)
efit <- eBayes(vfit)

summary(decideTests(efit,p.value = p.val_cutoff, lfc = log2(FC_cutoff)))

# These are the differential expression tables per cluster
C1 = topTable(efit, coef="C1_VS_ctrl", number = Inf)#,p.value = p.val_cutoff, lfc = log2(FC_cutoff))
A = topTable(efit, coef="A_VS_ctrl", number = Inf)#,p.value = p.val_cutoff, lfc = log2(FC_cutoff))

degList = list(A=A, C1 = C1)
saveRDS(degList, paste(diffExpr,".",dataMod,".subtypes.WINA.beta",beta,".RDS", sep = ""))
####################################### E N D ######################################