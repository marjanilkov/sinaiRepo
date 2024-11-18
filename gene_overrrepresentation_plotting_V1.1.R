rm(list = ls())

library(ggplot2)
library(tidyverse)
library(WGCNA)
library(SuperExactTest)
library(lattice)
source("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew3/20220829_data/codes_for_data/R-tomfunctions_static.R")


wd = ("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/20230816_miRNA_ad/v1/cor_div_2_av/freqThreshold5percent/")
setwd(wd)
cormethod = "pearson"
corrlpower = 1
#br.dn = read.delim("../csf_average_/CSF_brownd.tsv")
#turq.dn = read.delim("../csf_average_/CSF_turquoised.tsv")

n_path = 200 #Number of pathways
pval_cutoff = 0.05 # The cutoff value for the pvalue we want to use
# We will now create barplots with the pathway analyses

# Extract the sheet names which are the names of the experiments
library(readxl)
files_x0 = excel_sheets(path = "C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/20230816_miRNA_ad/v1/cor_div_2_av/freqThreshold5percent/enrichList.xlsx")
excel_data <- lapply(files_x0, read_excel, path = "C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/20230816_miRNA_ad/v1/cor_div_2_av/freqThreshold5percent/enrichList.xlsx")
names(excel_data) = files_x0

f_name = "red_pos_3"
filelist_x = list()
x = c() # The x column for the plotting
y = c() # the y column for the plotting
p.val = c() # The p.val column fot the plotting
ratio = c()
clusterGOdf = data.frame(matrix(NA,    # Create empty data frame
                                nrow = n_path,
                                ncol = 1))
colnames(clusterGOdf) = "X"
# the list of experiments needs to be fixed to take into account that some have NO SIGNIFICANT PATHWAYS
# files_x0 = c("red_pos_3","yellow_pos_1","yellow_neg_2","turquoise_neg_1","turquoise_pos_2"); projID = "sc1"
# files_x0 = c("red_neg_1","green_pos_3","blue_pos_4","yellow_neg_1","yellow_pos_2"); projID = "sc2"
# files_x0 = c("turquoise_pos_1", "turquoise_neg_2", "turquoise_pos_3"); projID = "sc3"
# files_x0 = c("blue_neg_3", "blue_pos_1", "blue_neg_4", "turquoise_neg_3"); projID = "sc4"
# files_x0 = c("red_pos_2", "red_neg_3"); projID = "sc5"

files_x = c()
for(f_name in files_x0)
{
  print(f_name)
  data1 = excel_data[[f_name]]
  tmp3 = data1[data1$P.adj < pval_cutoff,]
  if(nrow(tmp3)>0){
  #tmp3 = tmp3[1:n_path,]
  # we will change the clusters colors the same as in the paper A, B1, B2, C1, C2
  # We are extracting title for the plots from the file name
  z = gsub('tsv','',f_name)
  z = gsub('\\.',' ',z)
  z = gsub(" ", "", z)
  z = gsub("CSF_", "", z)
  
  tmp3$GOsets = gsub("^.{0,5}", "",  tmp3$GOsets) # to remove the GOBB or GOCB in the beginning of the sentence
  tmp3$GOsets = gsub("_", " ", tmp3$GOsets) # to remove all the underscores
  tmp3$GOsets = gsub("MORPHOGENESIS", "MORPH.", tmp3$GOsets)
  tmp3$GOsets = gsub("DEVELOPMENT", "DEV.", tmp3$GOsets)
  tmp3$GOsets = gsub("ACTIVITY", "ACT.", tmp3$GOsets)
  tmp3$GOsets = gsub("STRUCTURE", "STRUC.", tmp3$GOsets)
  tmp3$GOsets = gsub("ORGANIZATION", "ORG.", tmp3$GOsets)
  tmp3$GOsets = gsub("POSITIVE", "POS.", tmp3$GOsets)
  tmp3$GOsets = gsub("REGULATION", "REG.", tmp3$GOsets)
  tmp3$GOsets = gsub("EXTERNAL", "EXT.", tmp3$GOsets)
  tmp3$GOsets = gsub("BIOLOGICAL_ADHESION", "ADHESION.", tmp3$GOsets)
  tmp3$GOsets = gsub("EXTRACELLULAR", "EXTRACEL.", tmp3$GOsets)
  tmp3$GOsets = gsub("MULTICELLULAR", "MULTICEL.", tmp3$GOsets)
  tmp3$GOsets = gsub("EXTERNAL", "EXT.", tmp3$GOsets)
  tmp3$GOsets = gsub("SYSTEM", "SYS.", tmp3$GOsets)
  tmp3$GOsets = gsub("CIRCULATORY", "CIRC.", tmp3$GOsets)
  tmp3$GOsets = gsub("ANATOMICAL", "ANAT.", tmp3$GOsets)
  tmp3$GOsets = gsub("TRANSCRIPTION", "TRANSCR.", tmp3$GOsets)
  tmp3$GOsets = gsub("BIOSYNTHETIC", "BIOSYNT.", tmp3$GOsets)
  tmp3$GOsets = gsub("NEGATIVE", "NEG.", tmp3$GOsets)
  tmp3$GOsets = gsub("MEMBRANE", "MEMBR.", tmp3$GOsets)
  tmp3$GOsets = gsub("CONTAINING", "", tmp3$GOsets)
  tmp3$GOsets = gsub("SIGNALING", "SIGNAL.", tmp3$GOsets)
  tmp3$GOsets = gsub("MACROMOLECULE", "MACROMOL.", tmp3$GOsets)
  tmp3$GOsets = gsub("MODIFICATION_DEPENDENT", "", tmp3$GOsets)
  tmp3$GOsets = tolower(tmp3$GOsets)
  # truncate the pathways to only the first n letters
  #tmp3$GOsets = stringr::str_trunc(tmp3$GOsets, 40)
  
  filelist_x[[z]] = tmp3$GOsets
  x = c(x, rep((z), each=nrow(tmp3)))
  y = c(y, (tmp3$GOsets)[1:nrow(tmp3)])
  p.val = c(p.val, -log(tmp3$P.adj,10)[1:nrow(tmp3)])
  tmp3$ratio = (tmp3$GOsets.Size/tmp3$Overlap.Size)
  ratio = c(ratio, tmp3$ratio[1:nrow(tmp3)])
  #clusterGOdf = cbind(clusterGOdf, tmp3$GOsets)
  files_x = c(files_x, f_name)
  }
  }

goDF <- data.frame(x = x,
                   y = y,
                   p.val = p.val,
                   ratio = ratio)
# 
# colnames(clusterGOdf) = c("X", unique(goDF$x))
# clusterGOdf$X = NULL
# 
# 
# i =1; j = 1
# 
# tmp1 = matrix(NA, ncol = ncol(clusterGOdf), nrow = ncol(clusterGOdf))
# instead of trying to do average clustering we will use the overlap metric as a distance measurement
# for (i in 1:ncol(clusterGOdf))
# {
#   for(j in 1:ncol(clusterGOdf))
#   {
#     tmp1[i,j] = fisher.test(clusterGOdf[,i], clusterGOdf[,j])
#   }
# }

library(GOtest) ##from minghui
library(msigdb) ##from minghui
GOsets = c('c5.go.bp','c5.go.cc', 'c5.go.mf')
gosets_genes = msigdb.genesets(sets=GOsets, type='symbols', species='human',return.data.frame=T)
universe = curated.genesets(c('HGNC_universe'))$Gene
x = supertest(filelist_x, n = length(unique(gosets_genes$set)), degree = 2)

#openImgDev(paste(projID,"_msets.png", sep = ""),iwidth = 2048, iheight = 800, ipointsize = 32)
plot.msets(x, 
                      Layout = "landscape",
                      degree = 2, # the number of intersecting groups. It starts at 2 because otherwise it will show single gropus as intersections
                      keep.empty.intersections=FALSE,# whether to keep intersections with no overlap
                      sort.by=c('size','degree','p-value'))
#dev.off()
#xx = as.data.frame(x$set.names, x$P.value )
xx = summary(x)$Table

# Let's transform the pval array to a distance matrix
# We will only look at interactions of 1-on-1, not involving more or less than 2 
xx_2 = xx[lengths(regmatches(rownames(xx), gregexpr("1", rownames(xx))))==2,]
xx_2$cl1 = NA
xx_2$cl2 = NA
xx_2$cl1 = sub("\\&.*", "", xx_2$Intersections)
xx_2$cl2 = sub(".*\\&", "", xx_2$Intersections)
xx_2$cl1 = gsub(" ", "", xx_2$cl1)
xx_2$cl2 = gsub(" ", "", xx_2$cl2)

# Let's extract the cluster names
# we will change the clusters colors the same as in the paper A, B1, B2, C1, C2
# We are extracting title for the plots from the file name
z = gsub('tsv','',files_x)
z = gsub('\\.',' ',z)
z = gsub(" ", "", z)
z = gsub("CSF_", "", z)

dist_matrix = matrix(1, ncol = length(files_x), nrow = length(files_x), dimnames = list(z, z))
# ELt's populate the distance matrix with the values from the superExactTest
i = "turquoise_neg_1"
j = "turquoise_neg_2"

for (i in z)
{
  for (j in z)
  {
    if (i!=j)
    {
      print(i)
      print(j)
      # For distance we will use 1/(-log(pval,10))
      dist_matrix[i,j] = xx_2[(xx_2$cl1==j & xx_2$cl2==i) | (xx_2$cl2==j & xx_2$cl1==i),]$P.value
      print(".")
    }
    
  }
}

# To avoind 1/0 we will replace all the zeros with a small number
dist_matrix1 = ifelse(dist_matrix==1, 0.999, dist_matrix)
# We will now make the numbers in the distance matix be between 0 and 1 by 
# dividing everything by the largest number in the matrix
dist1 = 1+log(dist_matrix1,10)/max(-log(dist_matrix1,10))

h1row <- hclust(as.dist(dist1),method="complete")
collect_garbage()

# ----------------- output Hierarchical Clustering image -----------
#openImgDev("imgHierClust.png",iwidth = 2048, iheight = 800, ipointsize = 8)
par(mfrow=c(1,1))
plot(h1row, labels=F, xlab="",ylab="",main="",sub="")
#dev.off()

#myheightcutoff  ~ maximum dendrogam height allowed in modules
#mydeepSplit     ~ FALSE #TRUE will lead to iterative deep cut (resulting to many small modules), FALSE for normal dynamic cut
#myminModuleSize ~ minimal module size
mydeepSplit       = FALSE# fine structure within module
myminModuleSize   = 1 # modules must have this minimum number of genes
myheightcutoff = 0.95
mcolcode2= cutTreeStatic(hiercluster=h1row, heightcutoff=myheightcutoff, minsize1=myminModuleSize)

colcode.reduced  = reassignModuleNames(mcolcode2, minmodulesize=myminModuleSize, anameallmodules=F, 
                                       auseblackwhite=FALSE, useNumberAsLabel=FALSE, startlabel=1)
table(colcode.reduced)
plotDendrogramModuleLabels(mdendro=h1row, modulecolors=colcode.reduced, save2file=NULL, plotLabels=F)
mcex =  0.05 + 0.5/log10(ncol(dist1))
#openImgDev("imgCorHeatMap.png")
#openImgDev(paste(projID,"_heatmap.png", sep = ""),iwidth = 1500, iheight = 1024, ipointsize = 32)

#heatmapColorRG = ( rgcolors.func(50) )
heatmapColor = heat.colors(50)
par(mfrow=c(1,1))
diag(dist1) = 1.
heatmap(dist1,
        Rowv=as.dendrogram(h1row),
        Colv="Rowv", 
        scale="none",
        revC=T,
        #ColSideColors=as.character(colcode.reduced),
        #RowSideColors=as.character(colcode.reduced),
        cexRow = 1, 
        cexCol = 1,
        col=heatmapColor)
#dev.off()

################################################################################
# End clustering
################################################################################
# We will now arrange the columns of the data frame based on how they were clustered
clusterGOdf1 = data.frame(matrix(nrow = n_path, ncol = 1))

colnames(clusterGOdf1) = "X"
i = "turquoise_neg_1"
for  (i in names(filelist_x))
{
  tmp1 = filelist_x[[i]]
  tmp1 = tmp1[1:n_path]
  clusterGOdf1 = cbind(clusterGOdf1, tmp1)
}
colnames(clusterGOdf1)[2:ncol(clusterGOdf1)] = names(filelist_x)
clusterGOdf1$X = NULL

columns = data.frame(files_x, as.numeric(colcode.reduced))
columns = columns[order(columns$as.numeric.colcode.reduced.),]

new_goDF = goDF[0,]
clstName = "turquoise_pos_1"

for (clstName in columns$files_x)
{
  tmp1 = goDF[goDF$x == clstName, ]
  tmp1 = tmp1[1:n_path,]
  new_goDF = rbind(new_goDF, tmp1)
}

# theRearranger = rep(columns$colnames.clusterGOdf., each = n_path) 
# 
# goDF = goDF[match( goDF$x, theRearranger),]
# The probram arranges the x labels alphabetically, I am reordering them as they are
new_goDF$x <- factor(new_goDF$x, levels=unique(new_goDF$x))

ggplot(new_goDF, aes(x = new_goDF$x, y = new_goDF$y)) + 
  geom_point(aes(color = new_goDF$p.val, size = ratio)) #+
  #scale_size_continuous(range = c(5, 7)) # This line stops the smallest points 
#to become invisible instead it rescales everything to be between the values in range

# We will use ggplot2 to create the heatmap since ggplot2 offers more options 
# than the base plottong functions
library(reshape2)
#library(ggplot2)
meltDist = melt(dist1)
ggplot(meltDist, aes(Var1, Var2, fill= value)) + 
  geom_tile()+
  scale_fill_gradient(low="red", high="white") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Let's look at the frequency of which pathways appear in all the pathway sets
new_goDF1 = new_goDF[complete.cases(new_goDF),]
#library(dplyr)
#library(plyr)
new_goDF1 = new_goDF1%>% 
  add_count(y, name = "freq")
# We will create a new metric so as to try to ameliorate the following problem:
# a mildly significant pathway with p.adj~0.01 or so could be widely seen all over 
# while a very significant pathwaya with p.adj~1e-30 could appear in 70% of sets 
# and although biologically very meaningfull, it will be hidden. That is why we multiply 
# the frequency and the -log10(p.adj).
new_goDF1$metric = new_goDF1$freq*new_goDF1$p.val 
# Let's find out which pathways appear with the highes metric regardless of their set belonging
# new_goDF1$x = NULL
# new_goDF2 = ddply(new_goDF1,"y",numcolwise(sum))
new_goDF1 <- new_goDF1[order(new_goDF1$metric, decreasing=TRUE),]
plotDF <- new_goDF1[!duplicated(new_goDF1$y),][1:10,]

# add additional column which contains the angle by which each word will be 
# rotated in the wordcloud. As we can see here the angle is a multiple of 15 degrees
# from -15 to +15 degrees so that the words do not appear upside-down
plotDF <- plotDF %>%
  mutate(angle = 15 * sample(-1:1, n(), replace = TRUE, prob = c( 1, 1, 1)))
library(ggwordcloud)

#png(paste(projID,"_wordcloud.png", sep = ""),bg = "transparent")
ggplot(
  plotDF,
  aes(label = y, size = metric, color = metric, angle = angle), fill = "transparent") +
  geom_text_wordcloud(area_corr = F, rm_outside = TRUE) + # make the area of the whole word proportional to the value
  scale_size_area(max_size = 20) +
  theme_minimal() +
  theme(strip.text = element_text(size = 40, colour = "red")) +
  scale_color_gradient(high = "red", low = "purple")+
  guides(size = F)#+
  #facet_wrap(~type)
# dev.off()
# ggsave(
#   plot = p,
#   filename = paste(projID,"_wordcloud.png", sep = ""),
#   bg = "transparent",width = 50, height = 25, units = "cm"
# )
saveRDS(new_goDF1, paste(projID,".RDS", sep = ""))

