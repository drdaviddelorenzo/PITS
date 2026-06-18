## Script to create the example_cfr_data dataset for the PITS package.
## Run once from the package root: source("data-raw/create_data.R")
## Requires devtools / usethis.

example_cfr_data <- read.csv("data-raw/example_preintervention_data.csv",
                             stringsAsFactors = FALSE)
names(example_cfr_data) <- c("time", "outcome")

usethis::use_data(example_cfr_data, overwrite = TRUE)
