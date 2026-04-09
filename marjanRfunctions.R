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
# A function to find all the unique elements in each list object that appear in
# that particular object and in no other on the list. This one uses pure elements not like the one before

uniqListPURE = function(listObj){
  aunique = list()
  for (i in 1:length(listObj))
  {
    asub = (listObj[[i]]) # all the elements in the subtype
    tmp1 = listObj[-i]# all the unique elements in all 
    # the other subtypes
    for(j in names(tmp1))
    {
      tmp1[[j]] = (tmp1[[j]])
    }
    bsub = unique(unlist(tmp1)) 
    aunique[[paste(names(listObj)[i], ".uniq", sep = "")]] = setdiff(asub,bsub) # a set of unique genes that appear only in the 
    # first subtype and in no other 
  }
  return(aunique)
}
################################################################################
# The intersect() function finds the pairwise intersections between all vectors in a list
# and returns another list with the results
# multiIntersect = function(listObj){
#   intersectList = list()
#   for (i in names(listObj))
#   {
#     for (j in names(listObj))
#     {
#       if(i!=j)
#       {
#         intersectList[[paste(i,j, sep = ".")]] = intersect(listObj[[i]], listObj[[j]])
#       }
#     }
#   }
#   return(intersectList)
# }
multiIntersect = function(listObj){
  intersectList = list()
  for (i in 1:(length(listObj)-1))
  {
    for (j in (i+1):(length(listObj)))
    {
      if(i!=j)
      {
        intersectList[[paste(names(listObj[i]),names(listObj[j]), sep = ".")]] = intersect(listObj[[i]], listObj[[j]])
      }
    }
  }
  return(intersectList)
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
  
  tmp2.pval = tibble::deframe(z) # from dataframe to named vector
  sTest$P.value = tmp2.pval
  # expected overlap
  overlap.sizes = sTest$overlap.sizes[z$setOverlap]
  P.value = sTest$P.value[z$setOverlap]
  overlap.expected = sTest$overlap.expected[z$setOverlap]
  
  sTest1 = sTest
  
  sTest1$overlap.sizes = overlap.sizes
  sTest1$P.value = P.value
  sTest1$overlap.expected = overlap.expected
  }
  return(sTest1)
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
    tmp1.d = rownames(tmp1[tmp1$logFC < log2(1/FC.cutoff) & tmp1$adj.P.Val<p.val.cutoff,])
    deps.res[[paste(i, ".up", sep = "")]] = tmp1.u
    deps.res[[paste(i, ".dn", sep = "")]] = tmp1.d
  }
  return(deps.res)
}
# A function to take the DEGs split by contrast from a list and extract the names 
# of the DEGs according to p.cutoff and FC cutoff and directionality

depTableToNamesListNOadjPVAL = function(deps, FC.cutoff = 1.2, p.val.cutoff = 0.05)
{
  deps.res = list()
  for (i in names(deps))
  {
    tmp1 = deps[[i]]
    tmp1.u = rownames(tmp1[tmp1$logFC > log2(FC.cutoff) & tmp1$P.Value<p.val.cutoff,])
    tmp1.d = rownames(tmp1[tmp1$logFC < -log2(FC.cutoff) & tmp1$P.Value<p.val.cutoff,])
    deps.res[[paste(i, ".up", sep = "")]] = tmp1.u
    deps.res[[paste(i, ".dn", sep = "")]] = tmp1.d
  }
  return(deps.res)
}
# The same function as depTableToNamesList() but the DEGs are not split by directionality
depTableToNamesList.NO.DIR = function(deps, FC.cutoff = 1.2, p.val.cutoff = 0.05)
{
  deps.res = list()
  for (i in names(deps))
  {
    tmp1 = deps[[i]]
    tmp1 = rownames(tmp1[abs(tmp1$logFC) > log2(FC.cutoff) & tmp1$adj.P.Val<p.val.cutoff,])
    deps.res[[i]] = tmp1
  }
  return(deps.res)
}

# A function to change the names of DEGs in a list from A.dn, C1.up to 
# (t.A).DEG.dn and (t.C1).DEP.up etc
DEGnames = function(degs, dataType = "t", diffExpr = "DEG")
{
  for (i in 1:length(degs))
  {
    text = str_split(names(degs)[i], "\\.", n = Inf, simplify = FALSE)
    text[[1]][1] = paste("(",dataType,".", text[[1]][1], ").", sep = "")
    text[[1]][2] = paste(diffExpr, text[[1]][2], sep = ".")
    text = paste(text[[1]][1],text[[1]][2], sep = "")
    names(degs)[[i]] = text
  }
  return(degs)
}
# This is a way to fix overlapping legend title and numbers
# The legend had some overlaps between the text and the numbers so we are 
# moving the text a bit higher
# guides(fill = guide_colourbar(title.position = "top",
#                               title.vjust = 3,
#                               label.position = "left"))+



# A function that takes as na input the DEG tables in a list and produces a list 
# of degs based on FC and p.val, as well as their unique, and intersection genes
# This new version avoids commutations (intersect between a, b and b, a)
DEGmachine = function(deg, FC.cutoff = 1.2, p.val.cutoff = 0.05)
{
  ################################################################################
  # extract the unique symbols and intersections
  degList = depTableToNamesList(deg, FC.cutoff = FC.cutoff, p.val.cutoff = p.val.cutoff)

  # we need to keep the original deglist for later, so we make a copy
  degListTmp = degList
  # extract the unique genes and dump everithing in one list. When we search for unique
  # genes we only compare all the downreg gene sets among themselves and the 
  # upreg amongst themselves. We do not mix. Ask Bin why
  options(warn=-1)
  rm(uniq.d, uniq.u)
  options(warn=0)
  # We will not use any DEGs with contrasts involving MCI or NC to find unique symbols
  degListTmp = degListTmp[!grepl("MCI", names(degListTmp))]
  degListTmp = degListTmp[!grepl("NC", names(degListTmp))]

  # The way the situation is set up now we can have a contrast not just between 
  # cases and controls (ex. blue.ctrl) but also cases of one subtype VS cases of 
  # another subtype (ex. blue.turq). We don't want the latter case to be involved 
  # in the unique symbol identification because unique identification should give 
  # us the symbols unique to each subtype and that means DEGs extracted through 
  # only cases VS controls not case1 vs case2. So we will remove any situation 
  # where control is not a contrast
  degListTmp = degListTmp[grepl("ctrl", names(degListTmp), ignore.case = T)]
  # In this particular case we have only one such subtype.vs.subtype contrast and 
  # we will remove it directly
  #degListTmp = degListTmp[!grepl("blue.turq", names(degListTmp))]


  # Split the list into two lists of up and down regulated DEG/P sets
  uniq.d = degListTmp[grepl("\\.dn", names(degListTmp))]
  uniq.u = degListTmp[grepl("\\.up", names(degListTmp))]

  # extract the unique symbols  per direction
  uniq.d = uniqListPURE(uniq.d)
  uniq.u = uniqListPURE(uniq.u)
  names(uniq.d) = paste(names(uniq.d), ".dir", sep = "")
  names(uniq.u) = paste(names(uniq.u), ".dir", sep = "")
  # extract the unique symbols regardless of direction
  uniq = uniqListPURE(degListTmp)
  
  # We also want to have the intersections
  deg.intersect = multiIntersect(degListTmp)

  degList = c(degList, uniq, uniq.u, uniq.d, deg.intersect)
  # and finally, remove empty elements from the list
  degList = degList[lapply(degList,length)>0] ## you can use sapply,rapply
  # library(stringi)
  # x = stri_list2matrix(degList, byrow=F)
  # colnames(x) = names(degList)
  # openxlsx::write.xlsx(x, "t.DEP.xlsx")
  return(degList)
}

# an R package wrapper function to search PubMed for papers talking about 
# AD/Dementia/Ageing AND a given gene of interest. 
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
      pubQry <- get_pubmed_ids(paste(rowTerms[i],"[Text Word] AND ", colTerms[j],"[Text Word]", sep = ""))
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

################################################################################
# Sometimes ENSEMBL genes have a version at the end. When we have a list of 
# different sets of ENSEMBL genes with version attached, this function removes that
version.remove = function(list.object)
{
  for (i in names(list.object))
  {
    list.object[[i]] = gsub("\\..*","",list.object[[i]])
  }
  return(list.object)
}


deg.gene.2.ensembl = function(deg)
{
  for (i in names(deg))
  {
    tmp1 = deg[[i]]
    tmp1$symbol = rownames(tmp1)
    # Load a list to transform the ENSEMBL names to symbols and vice versa
    G_list = read.delim2("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/20240305_scz/data/knownGenes.geneid.tsv", header = F)
    colnames(G_list) = c("ENSEMBL.ID", "gene.symbol", "chr", "start", "end", "sign", "description")
    tmp1 = merge(tmp1, G_list, by.x = "symbol", by.y = "gene.symbol")
    rownames(tmp1) = tmp1$ENSEMBL.ID
    deg[[i]] = tmp1
  }
  return(deg)
}

# if there is a list with genes in ENSEMBL format it will transform them to Gene Symbol
list.ensembl.2.gene = function(gene.list)
{
  
  # ---- Progress bar: added lines ----
  n_items <- length(gene.list)
  pb <- utils::txtProgressBar(min = 0, max = n_items, style = 3)
  on.exit(close(pb), add = TRUE)  # ensure the bar closes even if an error occurs
  done <- 0L
  # ----------------------------------
  G_list = read.delim2("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/20240305_scz/data/knownGenes.geneid.tsv", header = F)
  colnames(G_list) = c("ENSEMBL.ID", "gene.symbol", "chr", "start", "end", "sign", "description")
  
  for ( i in names(gene.list))
  {
    tmp1 = as.data.frame(gene.list[[i]])
    colnames(tmp1)[1] = "ENSEMBLid"
    # Load a list to transform the ENSEMBL names to symbols and vice versa
    tmp1 = merge(tmp1, G_list, by.x = "ENSEMBLid", by.y = "ENSEMBL.ID")
    tmp1 = tmp1[!(is.na(tmp1$gene.symbol)),]
    gene.list[[i]] = tmp1$gene.symbol
    
    # ---- Progress bar: update counter ----
    done <- done + 1L
    utils::setTxtProgressBar(pb, done)
    # --------------------------------------
    
  }
  return(gene.list)
}

# Replace column names in df_data using df_map
replace_column_names <- function(df_data, df_map) {
  # Create a named vector for mapping
  ensembl_to_symbol <- setNames(df_map$Symbol, df_map$Geneid)
  
  # Replace column names if they exist in the mapping
  colnames(df_data) <- ifelse(
    colnames(df_data) %in% names(ensembl_to_symbol),
    ensembl_to_symbol[colnames(df_data)],
    colnames(df_data)  # keep original if no match
  )
  return(df_data)
}

# Pairwise enrichment between list elements (DE lists vs module lists)
pairwise_enrichment <- function(degs,
                                modules,
                                N = NULL,                 # numeric universe size
                                universe = NULL,          # character vector of background IDs
                                alternative = "greater",  # enrichment direction
                                adjust.method = "BH",
                                min_overlap = 1) {        # filter: only report pairs with >= this overlap
  
  # ---- Input checks ----
  if (is.null(N) && is.null(universe)) {
    stop("Provide either numeric N (background size) OR a universe vector of IDs.")
  }
  
  # Normalize names
  if (is.null(names(degs)))    names(degs)    <- paste0("degs_", seq_along(degs))
  if (is.null(names(modules))) names(modules) <- paste0("module_", seq_along(modules))
  
  # If a universe vector is supplied, clip sets to this background and derive N
  if (!is.null(universe)) {
    U <- unique(universe)
    N <- length(U)
    # Intersect sets with U to keep only valid background items
    degs    <- lapply(degs,    function(x) unique(intersect(x, U)))
    modules <- lapply(modules, function(x) unique(intersect(x, U)))
  } else {
    U <- NULL  # purely numeric background
    # Make sets unique to avoid double-counting overlaps
    degs    <- lapply(degs,    unique)
    modules <- lapply(modules, unique)
  }
  
  # ---- Computation ----
  out <- vector("list", length(degs) * length(modules))
  idx <- 1L
  
  for (i in seq_along(degs)) {
    A <- degs[[i]]
    a <- length(A) # size of DE set
    for (j in seq_along(modules)) {
      B <- modules[[j]]
      b <- length(B) # size of Module set
      k <- length(intersect(A, B)) # overlap
      
      # Skip low overlaps if requested
      if (k < min_overlap) {
        out[[idx]] <- data.frame(
          deg_set = names(degs)[i], module_set = names(modules)[j],
          a_deg = a, b_module = b, k_overlap = k, N_universe = N,
          expected_overlap = (a * b) / N,
          enrichment_ratio = ifelse((a * b) > 0, k / ((a * b) / N), NA_real_),
          fisher_p = NA_real_, hyper_p = NA_real_, odds_ratio = NA_real_
        )
        idx <- idx + 1L
        next
      }
      
      # 2x2 table components (must be >= 0)
      m11 <- k
      m12 <- a - k
      m21 <- b - k
      m22 <- N - a - b + k
      
      # Guard: if you only supplied numeric N, ensure the table is valid
      valid_table <- all(c(m11, m12, m21, m22) >= 0)
      
      if (!valid_table) {
        # This happens when N is too small relative to set sizes.
        # Prefer supplying a universe vector to avoid this.
        p_fisher <- NA_real_
        or <- NA_real_
      } else {
        mat <- matrix(c(m11, m12, m21, m22), nrow = 2, byrow = TRUE,
                      dimnames = list(c("DE", "NotDE"), c("Module", "NotModule")))
        ft <- suppressWarnings(fisher.test(mat, alternative = alternative))
        p_fisher <- ft$p.value
        or <- unname(ft$estimate)  # odds ratio
      }
      
      # Hypergeometric (equivalent to one-sided Fisher for enrichment)
      # P(X >= k) where X ~ Hypergeom(m=b successes in N, draw k=a)
      hyper_p <- if (valid_table) {
        phyper(q = k - 1, m = b, n = N - b, k = a, lower.tail = FALSE)
      } else {
        NA_real_
      }
      
      expected <- (a * b) / N
      out[[idx]] <- data.frame(
        deg_set = names(degs)[i],
        module_set = names(modules)[j],
        a_deg = a, b_module = b, k_overlap = k, N_universe = N,
        expected_overlap = expected,
        enrichment_ratio = ifelse(expected > 0, k / expected, NA_real_),
        fisher_p = p_fisher,
        hyper_p = hyper_p,
        odds_ratio = or
      )
      idx <- idx + 1L
    }
  }
  
  res <- do.call(rbind, out)
  # Multiple-testing correction across all pairs
  res$fisher_p_adj <- p.adjust(res$fisher_p, method = adjust.method)
  res$hyper_p_adj  <- p.adjust(res$hyper_p,  method = adjust.method)
  
  # Order by strongest enrichment (smallest Fisher p)
  res[order(res$fisher_p, na.last = TRUE), ]
}

# Copy paste this code to plot barplots
################################################################################
# Proportion of sex in subtypes + control
# p = ggplot(data=counts, aes(x=subtype, y=percentage, fill = sex, color = sex)) +
#   geom_bar(stat="identity", width = .7)+
#   scale_color_manual(values = c("black", "black"))+
#   scale_fill_manual(values=c("coral", "steelblue"))+
#   theme_classic()+
#   xlab("SCZ subtype")+
#   ylab("%")+
#   theme( axis.text.x = element_text(size=font.size, colour = "black"),
#          axis.text.y = element_text(size=font.size, colour = "black"),
#          text = element_text(size=font.size))
# ggsave('samples.percentage.png', 
#        egg::set_panel_size(p, 
#                            width=unit(2, "in"), 
#                            height=unit(2, "in")), 
#        width = 4, 
#        height = 4, 
#        units = 'in', 
#        dpi = 300)