rm(list = ls())

library(WGCNA)
library(reshape2)

wkdir = "C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew6/initial_cluster_out_test/";
setwd(wkdir)

ryanClassification = read.delim("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew2/20220816_data/Ryans_results_subtype_definitions_12_2019/meta.MSBB.subtypes.BM36.121919.tsv", sep = "\t")
ryanClassification = ryanClassification[,c("SynapseId", "ADsubtype", "ADsubtypeclass")]
protClassification = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew6/initial_cluster_out_test/prot.beta1/prot.2.subtypes.WINA.beta1.RDS")
methClassification = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew6/initial_cluster_out_test/meth.beta2/meth.subtypes.WINA.beta2.RDS")

# # For Ryan's cluster we also need to change for color coding the clusters to letters and numbers
ryanClassification$ADsubtype = ifelse(ryanClassification$ADsubtype=="turquoise" , "C1", ryanClassification$ADsubtype)
ryanClassification$ADsubtype = ifelse(ryanClassification$ADsubtype=="red" , "B1", ryanClassification$ADsubtype)
ryanClassification$ADsubtype = ifelse(ryanClassification$ADsubtype=="blue" , "B2", ryanClassification$ADsubtype)
ryanClassification$ADsubtype = ifelse(ryanClassification$ADsubtype=="orange" , "C2", ryanClassification$ADsubtype)
ryanClassification$ADsubtype = ifelse(ryanClassification$ADsubtype=="yellow" , "A", ryanClassification$ADsubtype)

# Remove controls and MCI
protClassification = protClassification[!(protClassification$cluster=="MCI" | protClassification$cluster=="control"),]
methClassification = methClassification[!(methClassification$cluster=="MCI" | methClassification$cluster=="control"),]
ryanClassification = ryanClassification[!(ryanClassification$ADsubtype=="MCI" | ryanClassification$ADsubtype=="control"),]

# For simplicity we will use the samples that appear in all three omics data types
# create a common ID to be able to properly knit together these DFs
commonID = Reduce(intersect,list(rownames(protClassification), ryanClassification$SynapseId,rownames(methClassification)))
protClassification = protClassification[commonID,]
methClassification = methClassification[commonID,]
ryanClassification = ryanClassification[ryanClassification$SynapseId %in% commonID,]

# Label the clusters by data source
ryanClassification$ADsubtype = paste("t.", ryanClassification$ADsubtype, sep = "")
protClassification$cluster = paste("p.", protClassification$cluster, sep = "")
methClassification$cluster = paste("m.", methClassification$cluster, sep = "")

colnames(protClassification)[2] = "prot.cluster"
colnames(methClassification)[2] = "meth.cluster"

outdf.p = merge(ryanClassification, protClassification, by.x = "SynapseId", by.y = "SynapseId")
outdf.m = merge(ryanClassification , methClassification, by.x = "SynapseId", by.y = "SynapseId")

rownames(outdf.p) = outdf.p$SynapseId
rownames(outdf.m) = outdf.m$SynapseId

# # Marjan's module's relabeled to best match Ryan's modules and then checked how
# # they overlap
modulesOverlap.prot = overlapTable(outdf.p$prot.cluster, outdf.p$ADsubtype)
modulesOverlap.meth = overlapTable(outdf.m$meth.cluster, outdf.m$ADsubtype)


# p.values
p.meth = modulesOverlap.meth$pTable
p.prot = modulesOverlap.prot$pTable
p.mashed = rbind(p.meth, p.prot)
pval = -log(p.mashed,10)
pval = round(pval, 2)
# Remove all pvalues larger than 0.05 in log form
pval = ifelse(pval<=-log(0.05,10), 0, pval)
############################################
# overlap number
o.meth = modulesOverlap.meth$countTable
o.prot = modulesOverlap.prot$countTable
o.mashed = rbind(o.meth, o.prot)

# Heatmap 
df.pval <- melt(pval) # Transform the matrix into long format
df.overlap <- melt(o.mashed) # Transform the matrix into long format
df = merge(df.pval, df.overlap, by = c("Var1", "Var2"))
# The axes tick need to contain the size of the sets. We will add this here
t.size = table(outdf.m$ADsubtype)
m.size = table(outdf.m$meth.cluster)
p.size = table(outdf.p$prot.cluster)

df$x.num = 0
df$y.num = 0

df$y.num = ifelse(df$Var2 == "t.A", t.size[1], df$y.num)
df$y.num = ifelse(df$Var2 == "t.B1", t.size[2], df$y.num)
df$y.num = ifelse(df$Var2 == "t.B2", t.size[3], df$y.num)
df$y.num = ifelse(df$Var2 == "t.C1", t.size[4], df$y.num)
df$y.num = ifelse(df$Var2 == "t.C2", t.size[5], df$y.num)

df$x.num = ifelse(df$Var1 == "m.A", m.size[1], df$x.num)
df$x.num = ifelse(df$Var1 == "m.C1", m.size[2], df$x.num)

df$x.num = ifelse(df$Var1 == "p.A", p.size[1], df$x.num)
df$x.num = ifelse(df$Var1 == "p.C1", p.size[2], df$x.num)


colnames(df) <- c("x", "y", "p.value", "overlap.size", "m.p.size", "t.size")
# df$x = gsub("\\.", "", df$x)
# df$y = gsub("\\.", "", df$y)
# manual order of axis ticks
# df$x <- factor(df$x, levels=c("m.A", "m.C1", "p.A", "p.C1"))
# df$y <- factor(df$y, levels=c("t.A", "t.B1", "t.B2", "t.C1", "t.C2"))
#library(ggh4x)

# Put the number of elements in parentheses next to the cluster and data source
df$x = paste(df$x, " (",df$m.p.size,")", sep = "")
df$y = paste(df$y, " (",df$t.size,")", sep = "")

ggplot(df, aes(x = x, y = y, fill = p.value)) +
  geom_tile(aes(fill = p.value)) + 
  geom_text(aes(label = overlap.size), size = 16) +
  coord_fixed()+
  scale_fill_gradientn(name = "-log10(p)", colours = c("white",rev(heat.colors(50))))+
  guides(fill = guide_colourbar(title = bquote(-log[10](p)), barwidth = 1.5, barheight = 15))+
  #geom_tile(color = "black", lwd = .5, linetype = 1)+
  theme(axis.text= element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 45),
        axis.title= element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 55), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  xlab("Other omics subtypes") + 
  ylab("Transcriptomic subtypes") +
  theme(text = element_text(size = 45)) #+
  # scale_x_discrete(guide = "axis_nested")+
  # scale_y_discrete(guide = "axis_nested")

ggsave(paste("prot.meth_overlap.ryan2.png", sep = ""), width = 30, height = 30, units = "cm")

