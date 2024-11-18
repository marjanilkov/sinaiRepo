
rm(list = ls())

library(dplyr)
library(stringr)


df<- read.csv("C:\\Users\\Marjan Ilkov\\OneDrive - The Mount Sinai Hospital\\Desktop\\MSSM\\20230817_csf_ad\\ADNI\\data\\old/CSF_noadj.data.csv")

# use only the baseline CSF data 
baseline<- df %>%
  arrange(PTID) %>%
  group_by(PTID) %>%
  slice(1) %>%
  ungroup()


############################################################################################################
## keep only the complete cases and run linear regression to regress out the effect of age, sex, and race
############################################################################################################


########## preprocessing ###################
# 1. remove empty or NA column in baseline 
your_data<- baseline
columns_to_check <- 17:ncol(your_data) 
filtered_data <- your_data[rowSums(!is.na(your_data[, columns_to_check]) & your_data[, columns_to_check] != "") == length(columns_to_check), ]

# 2. change character to numeric in remain column 
# Extract numeric values after ">" or "<" characters
columns_to_transform <- 17:ncol(filtered_data)
pattern <- paste0(c(">", "<"), collapse = "|")

transformed_values <- filtered_data %>%
  select(columns_to_transform) %>%
  mutate(across(.cols = everything(),
                .fns = ~ ifelse(grepl(pattern, .), str_extract(., "(?<=[><]).*"), .)))   


# Rename the columns to indicate the transformation
colnames(transformed_values) <- paste("transformed_", colnames(filtered_data)[columns_to_transform], sep = "")
filtered_data[colnames(transformed_values)]<- transformed_values
filtered_data[, 20:ncol(filtered_data)]<- apply(filtered_data[, 20:ncol(filtered_data)], 2, as.numeric)

########## adjustment ###################
### calculate adjusted value using linear regression
x<- filtered_data
i=20
for (i in 20:ncol(x)) {
  col_name <- colnames(x)[i]  # Get the column name
  AGE <- x$AGE
  GENDER <- x$Gender
  RACE<- x$RACE
  
  # Extract the current column and perform adjustment
  m_adjusted <- mean(x[[col_name]]) + residuals(lm(x[[col_name]] ~ AGE + GENDER + RACE))
  
  # Assign the adjusted column back to the dataframe
  x[[col_name]] <- m_adjusted
}

write.csv(x, 'CFS_adj.csv', row.names = F)