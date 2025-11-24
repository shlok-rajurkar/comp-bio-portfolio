# Functions to assist processing and cleaning dataset
required_packages <- c("dplyr", "tools")

for (package in required_packages){
  if (!(package %in% installed.packages()[, 1])) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}


available_fns <- c("store_file", "quantify_na_by_column",
                   "quantify_na_by_row", "remove_na_vals")
print("Available functions:")
for (fn in available_fns) {
  print(fn)
}

# Stores xpt or csv file as dataframe

store_file <- function(file_path) {
  read_fn <- {
    switch(tools::file_ext(file_path)
           , "xpt" = haven::read_xpt, "csv" = readr::read_csv)
  }
  dataframe <- read_fn(file_path)
  dataframe
}

## Following two functions can be used to quantify amount and distribution of
## missing values

# Returns a table with the ratio of each column that is NA

quantify_na_by_column <- function(dataframe) {

  if (!(inherits(dataframe, "data.frame") | inherits(dataframe, "tibble"))) {
    stop("Provide a dataframe.")
  }

  na_percents <- c()

  for (column in seq_len(ncol(dataframe))) {
    na_percent <- (sum(is.na(dataframe[[column]]))
                   / length(dataframe[[column]]))
    na_percents <- c(na_percents, na_percent)
  }

  column_names <- colnames(dataframe)

  na_frame <- data.frame("Column" = column_names, "NA Percents" = na_percents)

  na_frame
}

# Returns the ratio of rows that have no NA values

quantify_na_by_row <- function(dataframe) {
  if (!(inherits(dataframe, "data.frame") | inherits(dataframe, "tibble"))) {
    stop("Provide a dataframe.")
  }

  intact_count <- 0

  for (row in seq_len(nrow(dataframe))) {
    if (!(any(is.na(dataframe[row, ])))) {
      intact_count <- intact_count + 1
    }
  }
  intact_count / nrow(dataframe)
}

# Removes rows with any missing values and returns new dataframe or tibble

remove_na_vals <- function(dataframe) {
  if (!(inherits(dataframe, "data.frame") | inherits(dataframe, "tibble"))) {
    stop("Provide a dataframe.")
  }

  intact_rows <- dataframe %>% filter(complete.cases(.))

  intact_rows
}

# Functions to replace NA with given values
# like mean of column or random other value in column
# are readily available in R