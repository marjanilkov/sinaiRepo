
library(bayesbio)
library(easyPubMed)

extractPubMedIDs <- function(rowTerms, colTerms) {
## rowTerms: terms to be searched for rows of matrix, e.g. gene symbols, drug names 
## colTerms: terms to be searched for columns of matrix, e.g. disease names, conditions

l1 <- length(rowTerms)
l2 <- length(colTerms)
QryRes <- matrix("", l1,l2)
Counts <- matrix(0, l1,l2)
for (i in 1:l1) {
	for (j in 1:l2) {
		pubQry <- get_pubmed_ids(paste(rowTerms[i], colTerms[j]))
		QryRes[i,j] <- paste(unlist(pubQry$IdList),collapse="; ")
		Counts[i,j] <- pubQry$Count
	}
}
colnames(QryRes) <- paste(colTerms, "PMIDs", sep="_")
colnames(Counts) <- paste(colTerms, "Count", sep="_")
rownames(QryRes) <- rowTerms
rownames(Counts) <- rowTerms
return(data.frame(Counts,QryRes))
}
