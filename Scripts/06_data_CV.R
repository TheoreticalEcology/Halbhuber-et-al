###########################################################################
########################### Cross Validation ##############################
###########################################################################
library(readr)      
library(dplyr)
df <- readRDS("Data/df_complete_processed.rds")
df <- ungroup(df)

## I want to create 3 different Cross Validations 
# 1. random validation
# 2. exclude chemical names 
# 3. exclude species names


library(groupdata2)
set.seed(123)

# random cv
df$randomFold = sample.int(10, size = nrow(df), replace = T)

# chemical and species blocked
df$CAS = factor(df$CAS)
df$scientificNameStd = factor(df$scientificNameStd)

# chemical blocked 
chemicalFolds <- fold(
  data = df,
  k = 10,
  id_col = "CAS"
)

colnames(chemicalFolds)[231] = "chemicalFold"

# species blocked
speciesFolds <- fold(
  data = df,
  k = 10,
  id_col = "scientificNameStd"
)

colnames(speciesFolds)[231] = "speciesFold"

test = merge(df, chemicalFolds,sort = FALSE )
test1 = merge(test, speciesFolds,sort = FALSE ) # pro fold sind 175 (einmal 176) verschiedene CAS values 


# checken, reinschreiben oder sicherer: auch merge machen? 
colnames(test1)
#filtered_data <- test[test1$speciesFold == 5, ]
#sorted_rows<- sort(table(filtered_data$chemicalFolds))
#length(unique(filtered_data$scientificNameStd))

saveRDS(test1, file = "Data/df_final.rds")


### check class imbalance 
scientific_name_counts <- sort(table(df$scientificNameStd))
plot(scientific_name_counts)

scientific_name_proportions <- sort(prop.table(scientific_name_counts)*100)

data <- as.data.frame(sort(table(df_fish$scientificNameStd)))

data <- data[order(-data$Freq), ]
data$cumsum <- cumsum(data$Freq)
total_observations <- sum(data$Freq)
threshold <- 0.75 * total_observations
selected_species <- data[data$cumsum <= threshold, ]

pdf("Figures/speciesOf75.pdf", width = 17, height = 17)
par(mar = c(2,2,2,2))
par(cex = 2)
barplot(selected_species$Freq, 
        names.arg = selected_species$Var1,
        las = 2,
        xlim = c(0, 700), 
        main = "Species representing 75% of Observations", 
        xlab = "Frequency", 
        cex.names = 1,
        horiz = TRUE)

dev.off()


