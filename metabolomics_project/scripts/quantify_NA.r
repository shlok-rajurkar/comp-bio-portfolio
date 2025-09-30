# Takes in a dataframe and returns the percentage of NA values in each column

if (!("dplyr" %in% installed.packages()[, 1])) {
  install.packages('dplyr')
}

library("dplyr")

dataframe <- commandArgs(trailingOnly = TRUE)

if (class(dataframe) != 'data.frame') {
  stop('Provide a dataframe.')
}

