##################################################################
######## Assemly of dataset-chunks of EPA ECOTOX database ########
##################################################################
# source: https://cfpub.epa.gov/ecotox/search.cfm

# Dataframe for EPA ecotox search parameters 
search_parameters <- data.frame(
  Category = c("Habitat", "", "Chemicals", "", "Effect Measurements", "", "Endpoints", "", "Species", "", "Test Conditions", "", "Publication Options", ""),
  Parameter_Group = c("", "Aquatic", "", "", "Mortality Group", "", "Concentration Based Endpoints", "", "", "", "", "", "", ""),
  Name = c("", "Aquatic", "", "", "Mortality", "", "LC50", "", "", "", "", "", "Year Ending:", "Year Starting:"),
  Value = c("", "", "", "", "ALL", "", "ALL", "", "", "", "", "", "1970", "2024"),
  Additional_Info = c("", "", "", "", "", "", "", "", "", "", "", "", "", ""),
  stringsAsFactors = FALSE
)

# Adding the search run-time as a separate row
search_run_time <- data.frame(Category = "Search run-time", Parameter_Group = "", Name = "", Value = "2024-05-17 04:49:43", Additional_Info = "", stringsAsFactors = FALSE)
search_parameters <- rbind(search_parameters, search_run_time)


library(readr)
ECOTOX_data_1 <- read_delim("Data/EPA_ECOTOX/1915_1970.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)
ECOTOX_data_2 <- read_delim("Data/EPA_ECOTOX/1971_1974.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)
ECOTOX_data_3 <- read_delim("Data/EPA_ECOTOX/1975_1976.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)
ECOTOX_data_4 <- read_delim("Data/EPA_ECOTOX/1977_1979.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)
ECOTOX_data_5 <- read_delim("Data/EPA_ECOTOX/1980_1982.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)
ECOTOX_data_6 <- read_delim("Data/EPA_ECOTOX/1983_1985.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)
ECOTOX_data_7 <- read_delim("Data/EPA_ECOTOX/1987_1990.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)
ECOTOX_data_8 <- read_delim("Data/EPA_ECOTOX/1991_1992.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)
ECOTOX_data_9 <- read_delim("Data/EPA_ECOTOX/1993_1995.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)
ECOTOX_data_10 <- read_delim("Data/EPA_ECOTOX/1996_2000.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
ECOTOX_data_11 <- read_delim("Data/EPA_ECOTOX/2001_2008.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
ECOTOX_data_12 <- read_delim("Data/EPA_ECOTOX/2009_2017.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
ECOTOX_data_13 <- read_delim("Data/EPA_ECOTOX/2018_2024.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)

ECOTOX_data_all <- rbind(ECOTOX_data_1, ECOTOX_data_2, ECOTOX_data_3, ECOTOX_data_4, ECOTOX_data_5, 
                         ECOTOX_data_6, ECOTOX_data_7, ECOTOX_data_8, ECOTOX_data_9, ECOTOX_data_10, 
                         ECOTOX_data_11, ECOTOX_data_12, ECOTOX_data_13)

###
## Filter for 'CAS Number', 'Chemical Name', 'Species Scientific Name', 'Species Genus', 'Conc 1 Mean (Standardized)'  

ECOTOX_data_all_filtered <- ECOTOX_data_all[ , c(1,2,3,5, 6,12, 27),
                                             drop = FALSE]

#LC50 <- readRDS("Data/LC50_all.rds")


saveRDS(ECOTOX_data_all, file = "Data/LC50_all.rds")

saveRDS(ECOTOX_data_all_filtered, file = "Data/LC50_dataframe.rds")

