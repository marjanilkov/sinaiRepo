
rm(list = ls())
library(RCy3)
library(stringr)

FC.thresh = 1.5
p.val.thresh = 0.05

setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/OLD.20240305_scz_asd_ad/cm/scz/kda/")

# Load the Bayesian or MEGENA network
load("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/OLD.20240305_scz_asd_ad/cm/scz/kda/megena.res/MEGENA.Results.MSSM.PENN.PITT.cases.only.RData")
mn = el
rm(el)
mn = mn[,c("row", "col")]

# We will add a column to the mn that says whether the protein is up or down-regulated
clust_input <- readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/OLD.20240305_scz_asd_ad/cm/scz/wina.clust/initial_cluster_out_test_MPP/gene.beta1/DEG.gene.subtypes.WINA.beta1.MPP.RDS")
################################################################################
# The names of the genes in clust_input are in ENSEMBL IDs. We will replace them with gene symbols
# Replace ENSEMBL with the list from Vera
veraDF = read.delim2("knownGenes.geneid.tsv", header = F)
colnames(veraDF) = c("ENSEMBL", "symbol", "chr", "unknown1", "unknown", "unknown", "note")

i="blue"
for (i in names(clust_input))
{
  tmp1 = clust_input[[i]]
  tmp1$gene.name = gsub("\\..*","",rownames(tmp1))
  # Some genes (about 20 out of 19111 are duplicated. We will delete them)
  
  tmp1$gene.name =  veraDF$symbol[ match(tmp1$gene.name, veraDF$ENSEMBL) ]
  tmp1 = tmp1[!(duplicated(tmp1$gene.name)),]
  rownames(tmp1) = tmp1$gene.name 
  tmp1$gene.name = NULL
  clust_input[[i]] = tmp1
}
################################################################################
# Need to put the rownames as a column and add a directionality column and 
# filter for genes with p.val<0.05 and abs(FC)>FC,thresh
for (i in names(clust_input))
{
  tmp1 = clust_input[[i]]
  tmp1$gene_symbol = rownames(tmp1)
  tmp1$directionality = ifelse(tmp1$logFC>0, "up", "down")
  tmp1 = tmp1[abs(tmp1$logFC)>log2(FC.thresh) & tmp1$adj.P.Val<p.val.thresh,]
  clust_input[[i]] = tmp1
}

# We need to transform the data into SYMBOLS from ENSEMBL
# library('biomaRt')
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# G_list <- getBM(attributes= c("ensembl_gene_id", "hgnc_symbol", "uniprotswissprot"), mart= mart)
# # Some transripts have ENSEMBL symbol but for hgnc symbol they have an empty field , i.e. ""
# # We need to put their ENSEMBL symbol instead of the "" so as to not make errors afterwards
# G_list$hgnc_symbol = ifelse(G_list$hgnc_symbol == "", G_list$ensembl_gene_id, G_list$hgnc_symbol)

mn$row = gsub("\\..*","",mn$row)
mn$col = gsub("\\..*","",mn$col)

# Replace ENSEMBL with Glist from bioMart
# mn$row.symbol =  G_list$hgnc_symbol[ match(mn$row, G_list$ensembl_gene_id) ]
# mn$col.symbol =  G_list$hgnc_symbol[ match(mn$col, G_list$ensembl_gene_id) ]

mn$row.symbol =  veraDF$symbol[ match(mn$row, veraDF$ENSEMBL) ]
mn$col.symbol =  veraDF$symbol[ match(mn$col, veraDF$ENSEMBL) ]


# If an ensembl symbol has a gene symbol we put that, otherwise we keep the ensembl symbol
mn$row.symbol = ifelse(is.na(mn$row.symbol), mn$row, mn$row.symbol)
mn$col.symbol = ifelse(is.na(mn$col.symbol), mn$col, mn$col.symbol)

mn$row = NULL; mn$col = NULL
colnames(mn) = c("row", "col")
# Remove everything after the last |
# EXplanation of code
# s <- "DNS000001320_309.0/121.0_t0"
# t <- gsub("^([^_]*_[^_]*)_.*$", "\\1", s)
# t
# will print:
#   
#   DNS000001320_309.0/121.0
# A quick explanation of the regex:
#   
#   ^         # the start of the input
#   (         # start group 1
#     [^_]*   #   zero or more chars other than `_`
#       _       #   a literal `_`
#     [^_]*   #   zero or more chars other than `_`
#   )         # end group 1
# _         # a literal `_`
# .*        # consume the rest of the string
#   $         # the end of the input
#   which is replaced with:
#   
#   \\1       # whatever is matched in group 1
# And if there are less than 2 underscores, the string is not changed.
# mn$from = gsub("^([^\\|]*\\|[^\\|]*)\\|.*$", "\\1", mn$from)
# mn$to = gsub("^([^\\|]*\\|[^\\|]*)\\|.*$", "\\1", mn$to)
# # remove the sp or tr with :
# mn$from = gsub(":sp", "", mn$from)
# mn$to = gsub(":sp", "", mn$to)
# mn$from = gsub(":tr", "", mn$from)
# mn$to = gsub(":tr", "", mn$to)
# 
# # separate all words
# #Extract the middle word and use it as the symbol
# val = str_split(mn$from,"\\|")
# val = sapply(val,function(x) x[[2]])
# mn$source = val
# 
# val = str_split(mn$to,"\\|")
# val = sapply(val,function(x) x[[2]])
# mn$target = val

dataSet = mn[,c("row", "col")]  
sum(is.na(dataSet))
dataSet = dataSet[complete.cases(dataSet),]
colnames(dataSet) = c("source", "target")

# write.table(dataSet, "MPP.megena.v2.tsv", sep = "\t", row.names = F)
# createNetworkFromDataFrames(edges=dataSet)

# Let's find the subnetwork of the key drivers
#string.net<-getNetworkSuid()  #grab handle on network for future reference
KDAfile = read.delim("MEGENA-KeyDrivers-undirected-MPP.v2/WINA_subnets_AD_L1_KDx_combined.xls")
# replace the ensembl with gene symbol. But before that remove the version
KDAfile$keydrivers= gsub("\\..*","",KDAfile$keydrivers)
KDAfile$keydrivers1 =  veraDF$symbol[ match(KDAfile$keydrivers, veraDF$ENSEMBL) ]
# Since some ENSEBL symbols do not have gene symbols I will fill the NAs with the ensembl symbols again
KDAfile$keydrivers1 = ifelse(is.na(KDAfile$keydrivers1), KDAfile$keydrivers, KDAfile$keydrivers1)
KDAfile = KDAfile[complete.cases(KDAfile$keydrivers),]
KDAfile$keydrivers = KDAfile$keydrivers1
KDAfile$keydrivers1 = NULL
# Instead of splitting the subtype KDAs to up and down , we will combine them
subName = "turquoise"

#single.kdas = grep('SCZ', KDAfile$module, value=TRUE)
KDAs = KDAfile[KDAfile$module==subName,]$keydrivers
print(subName)

dirDF = clust_input[[subName]]
dirDF = dirDF[, c("gene_symbol", "directionality")]
#remove the version
dirDF$gene_symbol= gsub("\\..*","",dirDF$gene_symbol)
dirDF$symbol1 =  veraDF$symbol[ match(dirDF$gene_symbol, veraDF$ENSEMBL) ]
# Since some ENSEBL symbols do not have gene symbols I will fill the NAs with the ensembl symbols again
dirDF$symbol1 = ifelse(is.na(dirDF$symbol1), dirDF$gene_symbol, dirDF$symbol1)
dirDF$symbol = dirDF$symbol1
dirDF$gene_symbol = NULL
dirDF$symbol1 = NULL

dirDF = dirDF[!(dirDF$symbol == ""),]
dirDF = dirDF[complete.cases(dirDF$symbol),]
rownames(dirDF) = dirDF$symbol
# We will add a column stating whether a gene is a key driver or not
KDAs = as.data.frame(KDAs)
KDAs$is.kd = "T"
dirDF = merge(dirDF, KDAs, by.x = "symbol", by.y = "KDAs")

#write.csv(dirDF, paste("dirDF.scz.MPP.",subName,".casesonly.csv", sep = ""))

#-------------------------------------------------------------------------------
# We will transform the KDA names to the splice names. FInd all the possible places 
# where a symbol name appears, ex. for KRT 8 we get in FROM and TO "KRT8|F8W1U3" "KRT8|P05787"
# ast being the unique splice symbols appearing in the from and to columns
# i = "S100A9"
# KDAs.splice = c()
# for (i in KDAs)
# {
#   tmp1 = unique(mn$row[grepl( paste(i), mn$row, fixed = TRUE)])
#   tmp2 = unique(mn$col[grepl( paste(i), mn$col, fixed = TRUE)])
#   KDAs.splice = c(KDAs.splice, tmp1, tmp2)
# }

KDAs <- as.vector(as.matrix(KDAs))
KDAs = as.data.frame(unique(KDAs))

colnames(KDAs) = "KDAs"
# Now we will map these names to their short versions: Instead of "KRT8|F8W1U3" "KRT8|P05787"
# we will have KRT and KRT.1
# name.map = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/WGCNA_from_scratch/anew6/KDA/name.mapping.RDS")
# KDAs = merge(KDAs, name.map, by.x = "KDAs", by.y = "tmp_df")
# KDAs = KDAs$first_name2
#-------------------------------------------------------------------------------
KDAs = KDAs$KDAs
# The node.names needs to be a vector not a df
#mnSub = getFirstNeighbors(node.names = KDAs)#,network = string.net)
#get second neighbors
#bnSuba = getFirstNeighbors(node.names = bnSub)#,network = string.net)

# put together the KDAs and their first neightors for the network visualization
#mnSub1 = c(mnSub, KDAs)

selectNodes(KDAs, "name")
# GO to CYtoscape now and manually press CTRL+6 to select the next layer of nodes. 
# Each additional CTRL+6 selects one more layer of nodes
createSubnetwork(nodes = "selected")
# do this to create columns that map to node size and font sizes
analyzeNetwork(directed = T)

# To be able to show the whole MEGENA network I ned to label which genes belong 
# to which module
# replace the ensembl with gene symbol. But before that remove the version
module.output$id= gsub("\\..*","",module.output$id)
module.output$id =  veraDF$symbol[ match(module.output$id, veraDF$ENSEMBL) ]
sum(is.na(module.output$id))
write.csv(module.output, "module.output.csv")




#layoutNetwork('fruchterman-rheingold gravity_multiplier=1 nIterations=10')
#layoutNetwork("kamada-kawai defaultEdgeWeight=0.5 m_anticollisionSpringStrength=1000 m_averageIterationsPerNode=10000 m_nodeDistanceRestLengthConstant=1000")
# setNodeSizeMapping(
#   "Outdegree",
#   table.column.values = NULL,
#   sizes = c(1,100),
#   mapping.type = "c",
#   default.size = NULL,
#   style.name = NULL,
#   network = NULL)
# setNodeFontSizeMapping(
#   "Outdegree",
#   table.column.values = NULL,
#   sizes = NULL,
#   mapping.type = "c",
#   default.size = NULL,
#   style.name = NULL,
#   network = NULL)
#layoutNetwork("force-directed defaultSpringLength=700 defaultSpringCoefficient=3")

#getLayoutPropertyNames('kamada-kawai')
#getLayoutPropertyValue('kamada-kawai', "maxWeightCutoff")


#saveSession(paste(subName))
#exportNetwork(paste(subName))

#}

# testing. nothing of importance
b = dirDF
t = dirDF

b1 = b$symbol
t1 = t$symbol
b = KDAfile[KDAfile$module == "blue",]$keydrivers
t = KDAfile[KDAfile$module == "turquoise",]$keydrivers
sum(b1 %in% t1)
