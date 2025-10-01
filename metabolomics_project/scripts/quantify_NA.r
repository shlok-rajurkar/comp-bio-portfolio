# Takes in a dataframe and returns the percentage of NA values in each column

required_packages <- c("dplyr", "tools", "tidyverse")

for (package in required_packages){
  if (!(package %in% installed.packages()[, 1])) {
    install.packages(package)
  }
  library(package)
}

file_paths_entered <- commandArgs(trailingOnly = TRUE)

for (file_path in file_paths_entered) {
  quantify_NA(file_path)
}

quantify_na <- function(file_path) {
read_fn <- switch(file_ext(file_path), "xpt" = read_xpt, "csv" = read_csv)
dataframe <- read_fn(file_path)
if (!(class(dataframe) %in% c("data.frame", "tibble"))) {
  stop("Provide a dataframe.")
}

# na_means <- 

# for (column in ncol(dataframe) {

# })

#outMatrix <- matrix(data = dataVec, nrow = nrow, ncol = ncol, byrow = byrow)

}
