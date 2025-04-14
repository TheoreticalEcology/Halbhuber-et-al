# Download more than 10.000 datapoints from EPA ECOTOX database 
install.packages("ECOTOXr")
download_ecotox_data(
 target = "./"
)
#download_ecotox_data(tempdir())

library(ECOTOXr)
help(ECOTOXr)

list_ecotox_fields(
  which = c("default", "extended", "full", "all"),
  include_table = TRUE
)

search_ecotox(
  search,
  output_fields = list_ecotox_fields("default"),
  group_by_results = TRUE,
  compute = FALSE,
  as_data_frame = TRUE
)
