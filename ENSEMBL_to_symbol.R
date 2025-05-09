
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(attributes= c("ensembl_gene_id_version", "hgnc_symbol"), mart= mart)
saveRDS(G_list, "ensebleTOsymbol.List.RDS")
# If we need the ensembl version we can use 
 # G_list <- getBM(attributes= c("ensembl_gene_id_version", "hgnc_symbol"), mart= mart)


or 
library(org.Hs.eg.db)
annots <- select(org.Hs.eg.db, keys=mir129$`b_129-5p.p`, 
                 columns="SYMBOL", keytype="ENSEMBL")


# let's transfrom from UniProt to gene symbol for better visualization
my_protein_ids = unique(unlist(as.vector(mn)))
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # human

Uniprot = getBM(
  attributes=c('uniprotswissprot', 'hgnc_symbol'), 
  mart = mart)

colnames(Uniprot) <- c("Ensembl_ID", "UniProt" )
