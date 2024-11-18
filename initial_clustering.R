# ## R version 3.5.3-3.6.3
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")
# install.packages(c('RcppEigen', 'Rcpp'))
# install.packages(c("http://cran.nexr.com/src/contrib/Rclusterpp_0.2.3.tar.gz"), repos=NULL, type="source")
# BiocManager::install(c("impute", "pcaMethods"))
# install.packages(c('dynamicTreeCut', 'flashClust', 'foreach'))
#install.packages("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/tools/WINA-ryan-0.1.5.2.tar.gz", type="source",repos=NULL)

rm(list = ls())
### input params ###
source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/codes_for_data/R-tomfunctions_static.R")

library("WINA")
library("jsonlite")
library(WGCNA)
setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew6/")
wkdir = getwd()

num_cores = 7
beta = 1
clustering = "meth"
out_dir = paste(wkdir,"/initial_cluster_out_test/",clustering,".beta",beta,"/",sep="")
dir.create(out_dir, showWarnings = F)
setwd(out_dir)
#*-------------------------------------------------------------------------------------
#* STEP 0: read in gene information, expression data and consolidate data
#*
#*
#*-------------------------------------------------------------------------------------
# we have differing numbers of AD samples in the protein, methylation and RNAseq
# data, so we will do an intersection of samples with AD that appear 
# simultaneously in all data sets. We use the standardized (Z transformed) data
methDF = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/BM_36.MSBB_DMR270_mean_exp.CDR_age_sex_adjusted.std.RDS")
protDF = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/msbb.AMPAD.AD.proteomics.CDR_adjusted.uniq_genes.BM36.covcorr.ryan.std.RDS")
geneDF = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/BM_36.RNA.CDR_age_sex_adjusted.std.RDS")
meta = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/metadata.RDS")

if (clustering=="meth"){datExpr = as.data.frame(t(methDF))}
if (clustering=="prot"){datExpr = as.data.frame(t(protDF))}
if (clustering=="gene"){datExpr = as.data.frame(t(geneDF))}
rm(geneDF, protDF, methDF)
gc()

# Use only the AD patients
meta_AD = meta[meta$CDR>=1,]
datExpr = datExpr[rownames(datExpr) %in% meta_AD$SynapseBrainID,]
datExpr = as.matrix(datExpr)


RsquaredCut = 0.8
pcheicutoff = 0.5
minheight = 0.79
myminModuleSize = 1
print("oneplus=TRUE")
print("running WINA. here we GO...")
result_true = wina(datExpr,
              headCol=NA,
              outputDir=out_dir,
              fname='WINA',
              beta=beta,
              RsquaredCut=RsquaredCut,
              linkage="average",
              pcheicutoff=pcheicutoff,
              cormethod=c("pearson","spearman"),
              myminModuleSize=myminModuleSize,
              myheightcutoff=minheight,
              imagetype="png",
              gene.missing.rate.cutoff=0.5,
              sd.cutoff=NULL,
              sample.missing.rate.cutoff=0.5,
              impute=c('mean','knn'),
              compute.connectivity.statistics=TRUE,
              plot.heatmap=TRUE,
              heatmap.downsample.size=5000,
              heatmap.enhancefold=6,
              heatmap.useRaster=T,
              ncores=num_cores,
              oneplus=TRUE)



#plotDendrogramModuleLabels(mdendro=h1row, modulecolors=colcode.reduced, save2file=NULL, plotLabels=T)

print("saving Rdata...")
#save.image(file=paste(out_dir,"WINA_run.Rdata",sep=""))
print("saved.")
getwd()
# 
load("WINA_h1row_dendrogram.Rdata")
load("WINA_colcode-reduced_dendrogram.Rdata")
plot(h1row)
# Use Static TreeCut  instead of dynamic
#load("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew6/initial_cluster_out_test/meth.beta2/WINA_h1row_dendrogram.Rdata")
mydeepSplit       = FALSE# fine structure within module
myminModuleSize   = 3 # modules must have this minimum number of genes
myheightcutoff = 0.55
mcolcode2= cutTreeStatic(hiercluster=h1row, heightcutoff=myheightcutoff, minsize1=myminModuleSize)

colcode.reduced  = reassignModuleNames(mcolcode2, minmodulesize=myminModuleSize, anameallmodules=FALSE, 
                                       auseblackwhite=FALSE, useNumberAsLabel=FALSE, startlabel=1)
table(colcode.reduced)
plotDendrogramModuleLabels(mdendro=h1row, modulecolors=colcode.reduced, save2file=NULL, plotLabels=FALSE)
NewMatrix = as.data.frame(cbind(mcolcode2, as.character(colcode.reduced)))

NewMatrix$mcolcode2 = NULL

NewMatrix$SynapseId = rownames(NewMatrix)
colnames(NewMatrix)[1] = "meth.cluster"
NewMatrix = NewMatrix[,c("SynapseId", "meth.cluster")]

saveRDS(NewMatrix, paste("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew6/initial_cluster_out_test/meth.beta1/meth.subtypes.WINA.beta",beta,".RDS", sep = ""))
# we need to have the datExpr again to find which samples come from control and which from AD
methDF = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/BM_36.MSBB_DMR270_mean_exp.CDR_age_sex_adjusted.std.RDS")
protDF = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/msbb.AMPAD.AD.proteomics.CDR_adjusted.uniq_genes.BM36.covcorr.ryan.std.RDS")
geneDF = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/BM_36.RNA.CDR_age_sex_adjusted.std.RDS")
meta = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/metadata.RDS")

if (clustering=="meth"){datExpr = as.data.frame(t(methDF))}
if (clustering=="prot"){datExpr = as.data.frame(t(protDF))}
if (clustering=="gene"){datExpr = as.data.frame(t(geneDF))}
rm(geneDF, protDF, methDF)
gc()

classification = as.data.frame(cbind(h1row$labels, as.character(colcode.reduced)))
colnames(classification) = c("id", "cluster")
rownames(classification) = classification$id
# Now do visual inspection of the heatmap and dendrogram and do any merging if needed
#classification$cluster = ifelse(classification$cluster=="grey", "blue", classification$cluster)
#classification$cluster = ifelse(classification$cluster=="brown", "blue", classification$cluster)
tmp1 = datExpr[,1:2]
classification = merge(classification, tmp1, by.x="id", by.y="row.names", all.y = T)
# remove the two unnecessary columns
classification = classification[, 1:2]
rownames(classification) = classification$id

# Add the controls and the MCI
classificationRyan = read.delim("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew2/20220816_data/Ryans_results_subtype_definitions_12_2019/meta.MSBB.subtypes.BM36.121919.tsv", sep = "\t")
classification = merge(classificationRyan, classification, by.x = "SynapseId", by.y = "id", all.y = T)
classification$cluster = ifelse(classification$CDR == 0, "control", classification$cluster)
classification$cluster = ifelse(classification$CDR == 0.5, "MCI", classification$cluster)

# check if the clusters overlap with Ryan's original clusters and relabel them accordingly
overlapTable(classification$ADsubtype, classification$cluster)
classification$cluster = ifelse(classification$cluster=="blue", "A", classification$cluster)
classification$cluster = ifelse(classification$cluster=="turquoise", "C1", classification$cluster)
classification = classification[,c("SynapseId", "cluster")]
classification = classification[complete.cases(classification$cluster),]
rownames(classification)=classification$SynapseId
saveRDS(classification, paste("meth.subtypes.WINA.beta",beta,".RDS", sep = ""))

