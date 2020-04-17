#load packages
library(tidyverse)
library(purrr)

#load data - note: files are tab seperated and row.names must be set to null for
#protein_summary files
proteins_run1 <- read.table(file = '../proteomics/data/proteomics_run1_proteinsummary_proteinpilot_2019-10-25.txt', 
                            sep="\t", header=TRUE, row.names = NULL)

data_path <- "../proteomics/data/protein_pilot_export/"   # path to the data
files <- dir(data_path, pattern = "proteomics_run1") # get file names

# create a data frame holding the file names
data_list <- tibble(filename = files) %>% 
  mutate(file_contents = map(filename,                            # read files into
                             ~ read.table(file.path(data_path, .), sep="\t", header = TRUE, row.names = NULL)) # a new data column
  )  
data_list

# create single dataframe for filtering
df <- unnest(data_list) # give error about "cols"?

# check to see if reversed cols match and are therefore redundent
a <- na.omit(df$X115.114)
b <- na.omit(df$X114.115)
match(a, b)
#all values from the columns match and so only one is needed 

# select all of the ratio column names 
ratio_cols <- names(dplyr::select(df, matches('114|115|116|117')))
check <- str_subset(ratio_cols, "X")
check #12 combinations, should be 6

a <- str_subset(ratio_cols, "115.114")
b <- str_subset(ratio_cols, "116.114")
c <- str_subset(ratio_cols, "117.114")
d <- str_subset(ratio_cols, "117.116")
e <- str_subset(ratio_cols, "117.115")
f <- str_subset(ratio_cols, "116.115")

cols_to_remove <- c(a, b, c, d, e, f)

#attempt to filter out rows where these cols have data - not working 
test <- df %>% dplyr::filter(is.na(.[,cols_to_remove]))

# remove redundent cols 
df[, cols_to_remove] <- NULL

# check that the right number of cols are left 
ratio_cols <- names(dplyr::select(df, matches('114|115|116|117'))) #figure out how to remove rows where all of these cols are NA?
check <- str_subset(ratio_cols, "X")
check #6 combinations left 

# remove cols that only contain NAs
not_all_na <- function(x) any(!is.na(x))

df <- df %>%
  select_if(not_all_na)

num_df <- length(data_list$filename)

for (n in 1:num_df) {
  data_list$test[[n]] <- select_if(data_list$file_contents[[n]], not_all_na)
  
}

# remove X.Cov.95. cols with "RRRRRsp"  in them 
df$X.Cov.95. <- as.character(df$X.Cov.95.)

df_filtered <- df %>%
  dplyr::filter(!grepl("RRRRRsp", X.Cov.95.))

#return dataframe to list of dataframes to filter pvalues 
df_list <- df_filtered %>%
  nest(df = -filename)
# filter to rows where pvalue is less than 0.05


