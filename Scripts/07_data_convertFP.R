#################################################################################
# Data preparation for the tabular-data processing models 
# convert fingerprints as character strings to numeric representation as a matrix

df <- readRDS("Data/df_final.rds")

convert_to_vector <- function(fingerprint) {as.integer(unlist(strsplit(fingerprint, split = "")))}

fingerprints_list <- lapply(df$Fingerprint, convert_to_vector)

fingerprints_matrix <- do.call(rbind, fingerprints_list)

df <- cbind(df, fingerprints_matrix)


# give each fingerprint bit a categorical name (numeric not proccessable in ranger)
numbers <- as.character(1:166)
new_colnames <- paste0("bit", numbers)
colnames(df)[233:398] <- new_colnames

df <- as.data.frame(df)
df <- df %>% select (-Fingerprint, -CAS, -SMILES, -scientificNameStd, -meanConc)

saveRDS(df, file="Data/df_final&converted.rds")

# ## How many 0en are in each fold --> nearly 30 perscent
# sapply(1:10, function(fold) {
# mean(sapply(df_fish_rf[df_fish_rf$chemicalFold!=fold,-c(1:21, 902:905)], function(p) mean(p)) == 0)
# })
# table(df_fish_rf[df_fish_rf$chemicalFold!=1,-(1:21)]$bit831)
# mean(table(df_fish_rf[df_fish_rf$chemicalFold==1,-(1:21)]$bit831))
#
# ## How many 0en are in each fold --> nearly 30 perscent
# sapply(1:10, function(fold) {
#   mean(sapply(xc[df_fish_rf$chemicalFold!=fold,-c(1:21)], function(p) mean(p)) == 0)
# })
# table(df_fish_rf[df_fish_rf$chemicalFold!=1,-(1:21)]$bit831)
# mean(table(df_fish_rf[df_fish_rf$chemicalFold==1,-(1:21)]$bit831))
#
# hist(df_fish_rf[df_fish_rf$chemicalFold==1,5])
#
#
# dim(df_fish_rf[!duplicated(df_fish_rf$VEp),])