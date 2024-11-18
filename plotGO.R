# take the signatures from the drug repositioning project and do GO analysis and 
# try to see if there are any neuronal pathways

rm(list = ls())

library(GOtest) ##from minghui
library(msigdb) ##from minghui
library(ggplot2)

wkdir = "C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/20240109_emudra/v3/"
setwd(wkdir)

sigList = readRDS("sigListPreNat.2.40.RDS")

GOsets = c('c5.go.bp','c5.go.cc', 'c5.go.mf')
gosets_genes = msigdb.genesets(sets=GOsets, type='symbols', species='human',return.data.frame=T)
universe = curated.genesets(c('HGNC_universe'))$Gene

i = 3
# make it into a data frame for the same of the function working properly
tmp1 = sigList[[i]]$Symbol

# make it into a data frame for the same of the function working properly
tmp1 = data.frame(tmp1, "group")

colnames(tmp1) = c("gene","group")

result_weight1 = GOtest(x = tmp1,
                            go = gosets_genes,
                            query.population = unique(gosets_genes$genes),
                            background = 'query',
                            name.x = names(sigList)[i],
                            name.go = 'GOsets',
                            method = 'hypergeometric',
                            ncores = (detectCores(all.tests = FALSE, logical = TRUE) -1) )

write.table(result_weight1,paste(names(sigList)[i],".GO.tsv", sep = ""), sep="\t", quote=F, row.names=F)

# plot the GO terms
filelist_x = list.files(wkdir, pattern = ".tsv")
p.val_cutoff =0.05
n_pathways = 5
library(ggplot2)
library(scales)
library(forcats)

plotting = read.delim(f_name)
plotting = plotting[1:n_pathways,]
graph = ggplot(data = plotting, aes(x = fct_rev(fct_reorder(GOsets, -log10(P.adj))), y = -log10(P.adj), color = System)) + geom_col() + labs(x = NULL)
graph + coord_flip()

f_name = "preNat.GO.tsv"
# for clusters
#for(f_name in filelist_x)
#{
  print(f_name)
  data1 = read.delim(f_name)
  
  tmp2= data1[data1$Pvalue<p.val_cutoff,]
  tmp2 = tmp2[order(-tmp2$Overlap.Size),] # arrange in descending order
  tmp3 = tmp2[1:n_pathways,]
  # we will change the clusters colors the same as in the paper A, B1, B2, C1, C2
  # We are extracting title for the plots from the file name
  x = gsub('.tsv','',f_name)
  z = x
  #z = gsub('\\.',' ',x)
  #z = gsub(" ", " and ", z)
  #z = gsub("_", " ", z)
  #z = gsub("^d$", " down regulated", z)
  #z = gsub("^u$", " up regulated", z)
  
  tmp3$GOsets = gsub("^.{0,4}", "",  tmp3$GOsets) # to remove the GOBB or GOCB in the beginning of the sentence
  tmp3$GOsets = gsub("_", " ", tmp3$GOsets) # to remove all the underscores
  tmp3$GOsets = tolower(tmp3$GOsets)
  tmp3$GOsets = substr( tmp3$GOsets , start = 1 , stop = 40 )
  
  # ggplot(data=tmp3, aes(x=reorder(GOsets, Pvalue), y=Overlap.Size, fill = -log10(Pvalue))) +
  #   geom_bar(stat="identity")+
  #   ggtitle(z)+
  #   xlab("")+
  #   ylab("Overlap size")+
  #   scale_fill_distiller(palette = "YlOrRd", direction = 1)+
  #   scale_x_discrete(labels = label_wrap(30)) + # from the scales library
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  #         # axis.line=element_blank(),
  #         # #axis.text.x=element_blank(),
  #         # axis.text.y=element_blank(),
  #         # axis.ticks=element_blank(),
  #         # axis.title.x=element_blank(),
  #         # axis.title.y=element_blank(),
  #         panel.background=element_blank(),
  #         # panel.border=element_blank(),
  #         # panel.grid.major=element_blank(),
  #         # panel.grid.minor=element_blank(),
  #         plot.background=element_blank(),
  #         plot.margin=unit(c(10,10,10,10), "cm"),
  #         panel.spacing=unit(c(10,10,10,10), "cm"))
  # 
  
  # ggsave(paste(x,"X.png", sep = ""),
  #        width = 30, height = 30, units = "cm")
  ggplot(tmp3, aes(y=reorder(GOsets, -log10(Pvalue)), x=-log10(Pvalue), fill= "red")) +
    geom_col(position=position_dodge2(preserve='single'), color="black") +
    labs(x='-log10(p)', y='') +
    #scale_x_reverse() +
    #scale_y_discrete(labels=wrap_format(38))+
    theme_classic()+
    theme(legend.position = "none", text = element_text(size=40))#+
    #ggtitle("Prenat. to 2yr")
  ggsave(paste(x,"5paths.png", sep = ""), width = 30, height = 15, units = "cm")
################################################################################
# Plot the zscores
preNatScore = openxlsx::read.xlsx("preNat40.xlsx")
AiqunShortList <- openxlsx::read.xlsx("AiqunShortLIst.xlsx", colNames = F)
x = merge(AiqunShortList, preNatScore, by.x = "X1", by.y = "DrugName")
x = x[x$NormalizedScore>0,]
ggplot(x, aes(y=reorder(X1, NormalizedScore), x=NormalizedScore, fill= "red")) +
  geom_col(position=position_dodge2(preserve='single'), color="black") +
  labs(x='Normalized Score', y='') +
  #scale_x_reverse() +
  #scale_y_discrete(labels=wrap_format(38))+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=40))#+
#ggtitle("Prenat. to 2yr")
ggsave(paste("preNatDrugs40.png", sep = ""), width = 30, height = 15, units = "cm")

