#############################################
######### Finish data for modelling #########
#############################################
library(readr)
library(dplyr)

# In unserem Datensatz gibt es zu viele selben CAS Values für selbe Taxons aber mit unterschiedlichen Konzentrationen (toxizität) das ist verwirrend für das Modell, 
# deshalb berechnen wir den Mittelwert über diese Fälle und kommen dadurch zu einem wesentlich kleineren Datensatz 
# Ein Problem in dem entstehenden Datensatz ist, dass es einige CAS Values gibt mit nur einer Beobachtung, wir verwenden ihn aber trotzdem bis wir aus CAS eine kontinuierliche Variable gemacht haben 

df <- readRDS("Data/df_complete.rds")
# convert to continous variables 
convert_to_num <- function(df, columns) {
  for (col in columns) {
    df[[col]] <- as.numeric(df[[col]])
  }
  return(df)
}
columns <- c("BEl", "VEp", "REs", "OGp", "BLs", "PFs", "PFv", "CPt", "RMl") 
df <- convert_to_num(df, columns)

complete_dataset = df %>%  
  group_by(across(all_of(setdiff(colnames(df), "Concentration")))) %>% 
  summarise(var = var(Concentration), 
            max = max(Concentration), 
            meanConc = mean(Concentration),
  )

# Filters for number of CAS values > 1 
sorted_row_sums <- sort(table(complete_dataset$CAS))
un_int_sorted<- as.data.frame(sorted_row_sums)
top_ohne1 <- un_int_sorted[un_int_sorted[, "Freq"]>=2, ]
selected_cas_numbers_ohne1 <- as.vector(top_ohne1$Var1)
complete_dataset_ohne1 <- complete_dataset[complete_dataset$CAS %in% selected_cas_numbers_ohne1, ]
sort(table(as.numeric(complete_dataset_ohne1$CAS)))

nrow(unique(complete_dataset_ohne1[,c("CAS","scientificNameStd")])) 

# Set variable "var" to null 
complete_dataset_ohne1$var <- NULL
complete_dataset_ohne1$max <- NULL
#complete_dataset_ohne1_ohneNA <- na.omit(complete_dataset_ohne1)
# log transform to get a better response distribution 
complete_dataset_ohne1$logmeanConc = log(complete_dataset_ohne1$meanConc+1)
complete_dataset_ohne1 = complete_dataset_ohne1[complete.cases(complete_dataset_ohne1),]

# If NAs removed, new cases of single CAS values, hence exclude single CAS values again
sort(table(as.numeric(complete_dataset_ohne1$CAS)))
sorted_row_sums_ohne1_df <- sort(table(complete_dataset_ohne1$CAS))
un_int_sorted_ohne1_df <- as.data.frame(sorted_row_sums_ohne1_df)
top_ohne1_df <- un_int_sorted_ohne1_df[un_int_sorted_ohne1_df[, "Freq"]>=2, ]
selected_cas_numbers_ohne1_df <- as.vector(top_ohne1_df$Var1)
complete_dataset_ohne1_df <- complete_dataset_ohne1[complete_dataset_ohne1$CAS %in% selected_cas_numbers_ohne1_df, ]
sort(table(as.numeric(complete_dataset_ohne1_df$CAS)))
nrow(complete_dataset_ohne1_df)
#nrow(unique(complete_dataset_ohne1_df[,c("CAS","scientificNameStd")])) 


df1 = complete_dataset_ohne1_df
sort(table(df1$CAS))


saveRDS(df1, file = "Data/df_complete_processed.rds")


#############################################
########## only SMILES codes  ###############
#############################################
df <- readRDS("Data/SMILES_LC50_traits.rds")
convert_to_num <- function(df, columns) {
  for (col in columns) {
    df[[col]] <- as.numeric(df[[col]])
  }
  return(df)
}
columns <- c("BEl", "VEp", "REs", "OGp", "BLs", "PFs", "PFv", "CPt", "RMl") 
df <- convert_to_num(df, columns)

## In unserem Datensatz gibt es zu viele selben CAS Values für selbe Taxons aber mit unterschiedlichen Konzentrationen (toxizität) das ist verwirrend für das Modell, 
# deshalb berechnen wir den Mittelwert über diese Fälle und kommen dadurch zu einem wesentlich kleineren Datensatz (von ca. 5400 zu ca. 800 Beobachtungen)
# Ein Problem in dem entstehenden Datensatz ist, dass es einige CAS Values gibt mit nur einer Beobachtung, wir verwenden ihn aber trotzdem bis wir aus CAS eine kontinuierliche Variable gemacht haben 

complete_dataset_fish = df %>%  
  group_by(across(all_of(setdiff(colnames(df), "Concentration")))) %>% 
  summarise(var = var(Concentration), 
            max = max(Concentration), 
            meanConc = mean(Concentration),
            )

# Filters for number of CAS values > 1 
sorted_row_sums_ohne1 <- sort(table(complete_dataset_fish$CAS))
un_int_sorted_ohne1 <- as.data.frame(sorted_row_sums_ohne1)
top_ohne1 <- un_int_sorted_ohne1[un_int_sorted_ohne1[, "Freq"]>=2, ]
selected_cas_numbers_ohne1 <- as.vector(top_ohne1$Var1)
complete_dataset_ohne1 <- complete_dataset_fish[complete_dataset_fish$CAS %in% selected_cas_numbers_ohne1, ]
sort(table(as.numeric(complete_dataset_ohne1$CAS)))

# Set variable "var" to null 
complete_dataset_ohne1$var <- NULL
complete_dataset_ohne1$max <- NULL
#complete_dataset_ohne1_ohneNA <- na.omit(complete_dataset_ohne1)
# log transform to get a better response distribution 
complete_dataset_ohne1$logmeanConc = log(complete_dataset_ohne1$meanConc+1)
complete_dataset_ohne1 = complete_dataset_ohne1[complete.cases(complete_dataset_ohne1),]
# If NAs removed, new cases of single CAS values, thereby again, exclude single CAS values 
sort(table(as.numeric(complete_dataset_ohne1$CAS)))
sorted_row_sums_ohne1_df <- sort(table(complete_dataset_ohne1$CAS))
un_int_sorted_ohne1_df <- as.data.frame(sorted_row_sums_ohne1_df)
top_ohne1_df <- un_int_sorted_ohne1_df[un_int_sorted_ohne1_df[, "Freq"]>=2, ]
selected_cas_numbers_ohne1_df <- as.vector(top_ohne1_df$Var1)
complete_dataset_ohne1_df <- complete_dataset_ohne1[complete_dataset_ohne1$CAS %in% selected_cas_numbers_ohne1_df, ]
sort(table(as.numeric(complete_dataset_ohne1_df$CAS)))
nrow(complete_dataset_ohne1_df)

df1 = complete_dataset_ohne1_df
sort(table(df1$CAS))

saveRDS(df1, file = "Data/df_SMILES_processed.rds")

# Drop levels 
# df_eigenvectors_fish$CAS = droplevels(df_eigenvectors_fish$CAS)

