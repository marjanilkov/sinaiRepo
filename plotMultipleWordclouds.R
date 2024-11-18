# This is a script to look at pathways tsv files and find if any subtypes have
# overlapping pathways. it is a more of a syntactic analysis

rm(list = ls())

library(SuperExactTest)
library(wordcloud2) 
library(WGCNA)
library(ggplot2)
library(htmlwidgets)
library(ggwordcloud)

setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/papers/20230302_ad_multiomic_subtyping/fig6/")

filelist_x = list.files("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/papers/20230302_ad_multiomic_subtyping/fig6/", pattern = ".tsv")

# the number of pathways in the plots
n = 15

# List of all pathway analysis files divided by intersection (up-up regulated,
# down-down reg and up-down and down-up oppisitely regulated)
filelist_up = filelist_x[filelist_x %in% c('UPprot_Au.tsv', 
                                           'UPprot_C1u.tsv', 
                                           'UPprot_Xu.tsv')]

filelist_dn = filelist_x[filelist_x %in% c('DNprot_Ad.tsv', 
                                           'DNprot_C1.tsv', 
                                           'DNprot_Xd.tsv')]

filelist_op1 = filelist_x[filelist_x %in% c("prot_A.p.d.C1.p.u.tsv",
                                            "prot_A.p.u.C1.p.d.tsv")]
filelist_op2 = filelist_x[filelist_x %in% c("prot_A.p.d.X.p.d.tsv",
                                            "prot_A.p.u.X.p.u.tsv")]
filelist_op3 = filelist_x[filelist_x %in% c("prot_C1.p.d.X.p.u.tsv",
                                            "prot_C1.p.u.X.p.d.tsv")] 

i = 1
dataPath = data.frame()

# This for loop takes the set of files with appropriate names and prepares them 
# for the wordcloud
for (f_name in filelist_op3){
  print(f_name)
  df1 = read.delim(f_name)
  # Extract only the pathways that are significant taking into account adj.p.val 
  # not nominal p.val
  #data1 = data1[data1$P.adj<0.05,]
  # for the word cloud we will put -log(p,10) to replace freq
  df1 = df1[,c("GOsets", "P.adj")]
  df1$P.adj = ifelse(df1$P.adj==1, 0.99, df1$P.adj) # to remove the 0 values which 
  # create problem with division with zero later on
  df1$P.adj = -log(df1$P.adj,10)
  # add additional column which contains the angle by which each word will be 
  # rotated in the wordcloud. As we can see here the angle is a multiple of 15 degrees
  # from -15 to +15 degrees so that the words do not appear upside-down
  df1 <- df1 %>%
    mutate(angle = 15 * sample(-1:1, n(), replace = TRUE, prob = c( 1, 1, 1)))
  
  
  colnames(df1) = c("words", "freq", "angle")
  # We have to shorten some of the sentences
  df1$words = gsub("MORPHOGENESIS", "MORPH.", df1$words)
  df1$words = gsub("DEVELOPMENT", "DEV.", df1$words)
  df1$words = gsub("ACTIVITY", "ACT.", df1$words)
  df1$words = gsub("STRUCTURE", "STRUC.", df1$words)
  df1$words = gsub("ORGANIZATION", "ORG.", df1$words)
  df1$words = gsub("POSITIVE", "POS.", df1$words)
  df1$words = gsub("REGULATION", "REG.", df1$words)
  df1$words = gsub("EXTERNAL", "EXT.", df1$words)
  df1$words = gsub("BIOLOGICAL_ADHESION", "ADHESION.", df1$words)
  df1$words = gsub("EXTRACELLULAR", "EXTRACEL.", df1$words)
  df1$words = gsub("MULTICELLULAR", "MULTICEL.", df1$words)
  df1$words = gsub("EXTERNAL", "EXT.", df1$words)
  df1$words = gsub("SYSTEM", "SYS.", df1$words)
  df1$words = gsub("CIRCULATORY", "CIRC.", df1$words)
  df1$words = gsub("ANATOMICAL", "ANAT.", df1$words)
  df1$words = gsub("TRANSCRIPTION", "TRANSCR.", df1$words)
  df1$words = gsub("BIOSYNTHETIC", "BIOSYNT.", df1$words)
  df1$words = gsub("NEGATIVE", "NEG.", df1$words)
  df1$words = gsub("MEMBRANE", "MEMBR.", df1$words)
  df1$words = gsub("CONTAINING", "", df1$words)
  df1$words = gsub("SIGNALING", "SIGNAL.", df1$words)
  df1$words = gsub("MACROMOLECULE", "MACROMOL.", df1$words)
  df1$words = gsub("MODIFICATION_DEPENDENT", "", df1$words)
  #df1$words = gsub("NEURON_TO_NEURON", "NEURON2NEURON", df1$words)
  

  
  #
  df1$words = gsub("^.{0,5}", "", df1$words)
  df1$words = tolower(df1$words)
  df1$words <- gsub("_", " ", df1$words)
  # let's use on ly the first 50 pathways ordered in decreasing freq order
  df1 = df1[order(df1$freq, decreasing = T), ] 
  
  # We need to make the frequencies whole positive numbers
  # So we will divide them by the lowest value and extract the whole number
  df1$freq = as.integer(df1$freq)#/df1$freq[nrow(df1)])
 
  data1 = df1[1:n,]
  f_name = gsub(".tsv", "", f_name)
  data1$type = f_name
  
  dataPath = rbind(dataPath, data1)
  i = i+1
}
  # we will now add a legend
  # dataPath[nrow(dataPath) + 1,] <- list("a", max(dataPath$freq), 0,"legend")
  # dataPath[nrow(dataPath) + 1,] <- list("a", (max(dataPath$freq)+min(dataPath$freq))/2,0, "legend")
  # dataPath[nrow(dataPath) + 1,] <- list("a", min(dataPath$freq), 0,"legend")
  # 
  
  ggplot(
    dataPath,
    aes(label = words, size = freq, color = freq, angle = angle)) +
    geom_text_wordcloud(area_corr = F, rm_outside = TRUE) + # make the area of the whole word proportional to the value
    scale_size_area(max_size = 20) +
    theme_minimal() +
    theme(strip.text = element_text(size = 40, colour = "red"))+
    scale_color_gradient(high = "red", low = "purple")+
    guides(size = FALSE)+
    facet_wrap(~type)
  #ggsave("op_test.png", width = 90, height = 30, units = "cm")
  
  
#}
  
  

  