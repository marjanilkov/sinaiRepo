# A list of functions frequently used in R

# A function to find all the unique elements in each list object that appear in
# that particular object and in no other on the list

uniqList = function(listObj){
  aunique = list()
  for (i in 1:length(listObj))
  {
    asub = rownames(listObj[[i]]) # all the elements in the subtype
    tmp1 = listObj[-i]# all the unique elements in all 
    # the other subtypes
    for(j in names(tmp1))
    {
      tmp1[[j]] = rownames(tmp1[[j]])
    }
    bsub = unique(unlist(tmp1)) 
    aunique[[paste(names(listObj)[i], ".uniq", sep = "")]] = setdiff(asub,bsub) # a set of unique genes that appear only in the 
    # first subtype and in no other 
  }
  return(aunique)
}
################################################################################
# The intersect() function takes ONLY 2 objects in a list and returns their intersection
# this function uses the amazing function Reduce in R to find the intersect of 
# multiple sets-all of them at once 
multiIntersect = function(listObj){
  tmp1 = Reduce(intersect, listObj)
  return(tmp1)
}
################################################################################

# A function which takes as an input the results from the supertest() function 
# from the SuperExactTest library and produces a list of only the significant 
# overlaps containing names of elements in overlaps, number of elements and 
# elements in the overlap
library(stringr)
library(SuperExactTest)
sigIntersectList = function(sTest)
{
  # Extract the names of the sets significantly overlapping
  sTestSig = as.data.frame(sTest$P.value)
  colnames(sTestSig) = "P.value"
  sTestSig$barcode = rownames(sTestSig)
  sTestSig$set = deBarcode(sTestSig$barcode, sTest$set.names)
  sTestSig = sTestSig[complete.cases(sTestSig),]
  sTestSig = sTestSig[sTestSig$P.value<0.05,]
  
  i = 3
  j = 1
  
  intersectList = list()
  setNames = character()
  for (i in 1:nrow(sTestSig))
  {
    tmp1 = sTestSig$set[i]
    tmp2 = stringr::word(tmp1, 1:4, sep = "&")
    tmp2 = na.omit(tmp2)
    tmp2 <- gsub(" ", "", tmp2)
    
    # create a temp list that contains the requested sets in the significantly overlapping sets
    tmpList = list()
    for (j in 1:length(tmp2))
    {
      tmpList[[tmp2[j]]] = sTest$x[[tmp2[j]]]
    }
    tmp3 = Reduce(intersect, tmpList)
    intersectList[[tmp1]] = tmp3
  }
  return(intersectList)
}

################################################################################
# a function to make the names of GO terms more readable in a GO term list
# a function to make the pathway names shorter in GO term list
GOshort = function(files_x0)
{
  for(f_name in names(files_x0))
  {
    print(f_name)
    data1 = files_x0[[f_name]]
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
      files_x0[[f_name]] = tmp3
    }
  }
  return(files_x0)
}

# Intersect one OUTSIDE signature set with 10 subtype signatures
#plot only the ones involving DEG
# We will try to plot ony the mitoGenes combinations
SuperTestSign.VS.Subtype = function(sTest)
{
  overlap.sizes = sTest$overlap.sizes[c("00000000011",
                                        "00000000101",
                                        "00000001001",
                                        "00000010001",
                                        "00000100001",
                                        "00001000001",
                                        "00010000001",
                                        "00100000001",
                                        "01000000001",
                                        "10000000001")]
  P.value = sTest$P.value[c("00000000011",
                            "00000000101",
                            "00000001001",
                            "00000010001",
                            "00000100001",
                            "00001000001",
                            "00010000001",
                            "00100000001",
                            "01000000001",
                            "10000000001")]
  overlap.expected = sTest$overlap.expected[c("00000000011",
                                              "00000000101",
                                              "00000001001",
                                              "00000010001",
                                              "00000100001",
                                              "00001000001",
                                              "00010000001",
                                              "00100000001",
                                              "01000000001",
                                              "10000000001")]
  sTest$overlap.sizes = overlap.sizes
  sTest$P.value = P.value
  sTest$overlap.expected = overlap.expected
  return(sTest)
}

# Plot only the significant overlaps for the SuperExactTest
SuperTestOnlySignificant = function(sTest)
{
  # Make an exception when there are no significant values
  if (sum(sTest$P.value<0.05)>0)
  {
  tmp1 = sTest$P.value[sTest$P.value<0.05]
  z = t(as.data.frame(as.list(tmp1)))
  z = z[complete.cases(z),]
  z = as.data.frame(z)
  z$setOverlap = rownames(z)
  z$setOverlap = gsub("X", "",z$setOverlap)
  colnames(z)[1] = "p.val"
  z = z[,c("setOverlap", "p.val")]
  
  tmp1.pval = tibble::deframe(z) # from dataframe to named vector
  sTest$P.value = tmp1.pval
  # expected overlap
  overlap.sizes = sTest$overlap.sizes[z$setOverlap]
  P.value = sTest$P.value[z$setOverlap]
  overlap.Expected = sTest$overlap.expected[z$setOverlap]
  
  sTest1 = sTest
  
  sTest1$overlap.sizes = overlap.sizes
  sTest1$P.value = P.value
  sTest1$overlap.Expected = overlap.Expected
  }
  return(sTest)
}

# A function to take the DEGs split by contrast from a list and extract the names 
# of the DEGs according to p.cutoff and FC cutoff and directionality

depTableToNamesList = function(deps, FC.cutoff = 1.2, p.val.cutoff = 0.05)
{
  deps.res = list()
  for (i in names(deps))
  {
    tmp1 = deps[[i]]
    tmp1.u = rownames(tmp1[tmp1$logFC > log2(FC.cutoff) & tmp1$adj.P.Val<p.val.cutoff,])
    tmp1.d = rownames(tmp1[tmp1$logFC < -log2(FC.cutoff) & tmp1$adj.P.Val<p.val.cutoff,])
    deps.res[[paste(i, ".UP", sep = "")]] = tmp1.u
    deps.res[[paste(i, ".DOWN", sep = "")]] = tmp1.d
  }
  return(deps.res)
}