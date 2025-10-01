# Takes in a dataframe and returns the percentage of NA values in each column

required_packages <- c("dplyr", "tools")

for (package in required_packages){
  if (!(package %in% installed.packages()[, 1])) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

file_paths_entered <- commandArgs(trailingOnly = TRUE)

for (file_path in file_paths_entered) {
  print(file_path)
  print(quantify_NA(file_path))
}

quantify_na <- function(file_path) {
  read_fn <- {
    switch(tools::file_ext(file_path)
           , "xpt" = haven::read_xpt, "csv" = readr::read_csv)
  }
  dataframe <- read_fn(file_path)
  if (!(inherits(dataframe, "data.frame") | inherits(dataframe, "tibble"))) {
    stop("Provide a dataframe.")
  }

  na_percents <- c()

  for (column in seq_len(ncol(dataframe))) {
    na_percent <- (sum(is.na(dataframe[[column]])) / length(dataframe[[column]])) * 100
    na_percents <- c(na_percents, na_percent)
  }

  column_names <- colnames(dataframe)

  na_frame <- data.frame("Column" = column_names, "NA Percents" = na_percents)

  na_frame
}
