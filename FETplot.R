# Create Fisher's exact test plot
rm(list = ls())

setwd("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/20240305_scz_asd_ad/cm/scz/centerEnrichmentCheck/")

meta = readRDS("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/2024/20240305_scz_asd_ad/cm/scz/Release4/Metadata/METADATA.with.CLUSTERS.RDS")

# lets load the classification
#remove the controls, NCs, and  MCI
meta2 = meta[meta$Dx == "SCZ",]

# Make 1000 fake genes
samples <- meta2$SampleID

# Make a bunch of fake DEG data
# DEGs for two subtypes
subtype_samples <- meta2[,c("SampleID", "cluster")]
# DEGs from 3 cohorts
cohort_samples <- meta2[,c("SampleID", "Institution")] # Cohort labels

#### Perform hypergeometric test ####
hyper <- GOtest(x = subtype_samples, 
                go = cohort_samples, 
                query.population = samples, 
                background = "query",
                name.x = "cluster", 
                name.go = "Institution", 
                method = "hypergeometric")
#


















phyper(q = 49 - 1, m = 153, n = 97, k = 141, lower.tail = F, log.p = F)


sum(meta2$Institution=="Pitt" & meta2$cluster == "blue")
phyper(q = 27 - 1, m = 97, n = 153, k = 141, lower.tail = F, log.p = F)

sum(meta2$Institution=="MSSM" & meta2$cluster == "blue")
phyper(q = 49 - 1, m = 97, n = 153, k = 141, lower.tail = F, log.p = F)

# Let's compare MSSM vs Pitt because MSSM are heavy cases which were 
# institutionalized followed for a long time while Pitt are post mortem which
# had SCZ without following and most are not institutionalized
meta3 = meta2[meta2$Institution == "MSSM" | meta2$Institution == "Penn",]
# create a contingency table
dat = table(meta2[,c("Institution", "cluster")])

# draw a mosaic plot for visualization
mosaicplot(dat, main = "", color = TRUE)
ggsave("FET.Institution.png", width = 7.5, height = 5)

chisq.test(dat)
test = fisher.test(dat)
fisher.test(table(meta3$cluster, meta3$Institution))
# Let's make a nicer plot
# combine plot and statistical test with ggbarstats
library(ggstatsplot)
ggbarstats(
  meta3, Institution, cluster,
  results.subtitle = F, 
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    round(test$p.value, 4)
  )
) + scale_fill_manual(values = c("white", "grey"))
#ggsave("FET.Institution.MSSM.Penn.png", width = 7.5, height = 5)


