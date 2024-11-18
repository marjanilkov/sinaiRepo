# a short script to compare the old and new ROSMAP subtypes and classes

rm(list = ls())

old = read.delim("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/ROSMAP.subtypes.matchedMSBB.final.txt")
new = read.delim("C:/Users/Marjan Ilkov/OneDrive - The Mount Sinai Hospital/Desktop/MSSM/ROSMAP_expanded_new_subtype_predictions_042423_v2.csv", sep = ",")

new = new[,c("projid", "large_types", "cogdx2")]
all = merge(old, new, by = "projid", all = T)
