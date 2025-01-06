# A script to produce Venn diagrams of DEGs
rm(list = ls())

setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/papers/20230302_ad_multiomic_subtyping/fig3")
source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/sinaiRepo/marjanRfunctions.R")

library(ggVennDiagram)
library(ggrepel)

################################################################################
# User defined vars
FC_cutoff = 1.2
p.val_cutoff = 0.1
degList = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/papers/20230302_ad_multiomic_subtyping/fig3/DEP.t.sub.RYAN.RDS")
## Use either depTableToNamesList or depTableToNamesListNOadjPVAL depending 
# whether the pval threshold is for the corrected or the nominal pval
degList = depTableToNamesListNOadjPVAL(degList, FC.cutoff = FC_cutoff, p.val.cutoff = p.val_cutoff)
# DOWN or UP regulated DEGs?
tmpChar = "dn"
dataType = "t"
dataMod = "DEP"
set.Names = "n" # whether to include the set names or not. options are y/n
################################################################################

tmpList = degList[grepl( tmpChar, names(degList), fixed = TRUE)]
degs = DEGnames(degs = tmpList, dataType = dataType, diffExpr = dataMod)

if (set.Names == "y"){NN = names(degs)} else{NN = ""}

ggVennDiagram(degs,  
              category.names = NN,
              label_alpha = 0,
              label = "count", 
              label_size = 15,
              set_size = 30) + 
  theme(text = element_text(size = 30),
        legend.title = element_text(size=50),
        legend.key.height= unit(2, 'cm'),
        legend.key.width= unit(2, 'cm'),
        legend.position = "none",
        legend.text = element_text(colour="black", size=40),
        plot.title = element_text(size = 60, face = "bold"))+
  scale_fill_gradient(low = "#F4FAFE", high = ifelse(tmpChar == "dn",  "#4981BF", "red"))+
  scale_x_continuous(expand = expansion(mult = .23))

#ggsave(paste(dataMod,tmpChar,".Venn.",ifelse(set.Names=="y", "YESnames", "NOnames" ),".png", sep = ""), width = 50, height = 40, units = "cm",limitsize = F)

