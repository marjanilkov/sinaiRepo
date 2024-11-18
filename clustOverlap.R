rm(list = ls())

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
options(stringsAsFactors=F)
library(DESeq2)
library(limma)
library(edgeR)
library(qvalue)
library(SuperExactTest)
set.seed(12345)
library(GOtest) ##from minghui
library(msigdb) ##from minghui
library(tidyverse)

source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/marjanRfunctions.R")
################################################################################
# User defined vars
degList = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew6/DEG/DEG.transcriptomic.subtypes.RYAN.RDS")
wkdir = "C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/papers/20230302_ad_multiomic_subtyping/fig2/";
p.val_cutoff =0.05
dataMod = "DEG"
modSize = nrow(degList[[1]])
dataType = "t"
FC = 1.2
# Use either depTableToNamesList or depTableToNamesListNOadjPVAL depending 
# whether the pval threshold is for the corrected or the nominal pval
degs = depTableToNamesList(degList,  p.val.cutoff = p.val_cutoff, FC.cutoff = FC)
################################################################################
setwd(wkdir)

# WE will rearrange the order of the DEG sets for nicer view of the axes
#degs = degs[c("A.dn", "B1.dn", "B2.dn", "C1.dn", "C2.dn", "A.up", "B1.up", "B2.up", "C1.up", "C2.up")]

degs = DEGnames(degs = degs, dataType = dataType, diffExpr = dataMod)
sTest = supertest(degs, n = modSize, degree = 2)

png(file = paste(dataMod,".png", sep = ""),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6, units = "in", res = 300) # The height of the plot in inches
plot.msets(sTest,
           Layout = "landscape",
           degree = 2:10, # the number of intersecting groups. It starts at 2 because otherwise it will show single gropus as intersections
           keep.empty.intersections=F,# whether to keep intersections with no overlap
           sort.by=c('size','degree','p-value'))#, 
# cex=2, cex.axis = 2, legend.text.cex = 2, overlap.size.cex = 2, 
# color.scale.cex = 1.5, color.scale.pos = c(0.75, 0.75))

dev.off()
# OVERLAP HEATMAP
library(GeneOverlap)
data(GeneOverlap)

gom.obj <- newGOM(degs, degs, modSize)

x = getMatrix(gom.obj, name="pval")
diag(x) = 1

x = -log(x,10)
x[x==Inf] <- 300
diag(x) = 0

index_list = which(x<0.05 , arr.ind = TRUE)
nList = getNestedList(gom.obj, name = "intersection") # list of genes inside intersections

################################################################################
# Plotting the heatmap of pvalues

library(reshape2)
x = melt(x) 
# To increase the visibility of the colors we will set the max at 20
#x$value = ifelse(x$value>20, 20, x$value)

ggplot(x, aes(Var1, Var2, fill= value)) + 
  geom_tile()+
  scale_fill_gradientn(name = "-log10(p)",
                        colours = c("white",rev(heat.colors(50))))+ # the 0 is
  # white and the rest of the colors are the same we use in the adjacency matrix 
  # but reversed
  geom_tile(color = "black",
            lwd = .5,
            linetype = 1) +
  geom_segment(aes(x = 5.5, y = 0.4, xend = 5.5, yend = 10.6), size = 6, linetype = 1)+
  geom_segment(aes(x = 0.4, y = 5.5, xend = 10.6, yend = 5.5), size = 6, linetype = 1)+
  
  # We draw two rectangles to outline the up and down regulated intersections
  # geom_rect(aes(xmin = 0.5, xmax = 5.5, ymin = 0.5, ymax = 5.5), 
  #           fill = "blue", alpha = 0., color = "blue", linewidth = 2)+
  # geom_rect(aes(xmin = 5.5, xmax = 10.5, ymin = 5.5, ymax = 10.5), 
  #           fill = "red", alpha = 0., color = "red", linewidth = 2)+
  guides(fill = guide_colourbar(title.position = "top",
                                title.vjust = 3,
                                label.position = "left"))+
  # theme_void()+
  coord_fixed()+
  theme(legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=40), #change legend title font size
        legend.text = element_text(size=30),
        plot.title = element_text(size=50),
        axis.text.x = element_text(face="bold", color="black", 
                                   size=34, angle = 45,
                                   hjust = 1,
                                   vjust = 1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=34))
  #ggtitle(titleText)

ggsave(paste(dataMod,"_clustOverlap.",p.val_cutoff,".png", sep = ""), width = 40, height = 30, units = "cm")

