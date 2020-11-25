# Note: this script is redundent, has be replaced by function

#load packages
library(tidyverse)
library(purrr)

#load data - note: files are tab seperated and row.names must be set to null for
#protein_summary files

# note: data from proteinpilot software contains empty cols and is a tab
# seperated .txt file - if those are loaded in as they are the data is not
# formated corretely - fixed by openeing in spreadsheet software, deleting empty
# cols and saving as a .csv
data_path <- "../proteomics/data/protein_pilot_export/empty_cols_removed/"   # path to the data
files <- dir(data_path, pattern = "proteomics_run2") # get file names

# create a data frame holding the file names
data_list <- tibble(filename = files) %>% 
  mutate(file_contents = map(filename,                            # read files into
                             ~ read.table(file.path(data_path, .), sep="\t", header = TRUE, row.names = NULL)) # a new data column
  )  
data_list


# remove cols that only contain NAs
#function below returns TRUE if there are only NAs in an object
not_all_na <- function(x) any(!is.na(x))


num_df <- length(data_list$filename)
# below uses select_if from dplyr to only select cols that only contain NAs
for (n in 1:num_df) {
  data_list$df_filtered[[n]] <- select_if(data_list$file_contents[[n]], not_all_na)
  
}
data_list
#convert column to character for filtering 
str(data_list$df_filtered[[1]])

for (n in 1:num_df) {
  data_list$df_filtered[[n]]$Accession <- as.character(data_list$df_filtered[[n]]$Accession)
}



str(data_list$df_filtered[[1]])


# remove X.Cov.95. cols with "RRRRRsp"  in them 
data_list
for (n in 1:num_df) {
  data_list$df_filtered[[n]] <- data_list$df_filtered[[n]] %>%
    dplyr::filter(!grepl("RRRRRsp", Accession))
}
data_list

# pval_cols <- list()
# for (n in 1:num_df) {
#   pval_cols[[n]] <- names(dplyr::select(data_list$df_filtered[[n]], contains("PVal")))
# }
# pval_cols
# 
# # filter to rows where pvalue is less than 0.05 - not working, dataframes have no observations
# test <- list()
# for (n in 1:num_df) {
#   test[[n]] <- data_list$df_filtered[[n]] %>%
#     dplyr::filter_at(vars(contains("PVal")), all_vars(.< 0.05))
# }
# str(test[[1]])



# create single dataframe for filtering
df_list <- data_list
df_list$file_contents <- NULL
#remove extra file 
#df_list <- df_list[-c(5),] 
df_list

# this unnest is not working correctly 
df <- unnest(df_list)

# check to see if reversed cols match and are therefore redundent - note: this
# used to be true when loading in the incorrect .txt file for reasons I don't
# understand - a ratio won't give the same value if you reverse it to this never
# should have worked
#a <- na.omit(df$X115.114)
#b <- na.omit(df$X114.115)
# below will return TRUE if any value doesn't match
#any(!match(a, b))
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
cols_to_remove


# remove redundent cols 
df[, cols_to_remove] <- NULL

# check that the right number of cols are left 
ratio_cols <- names(dplyr::select(df, matches('114|115|116|117'))) #figure out how to remove rows where all of these cols are NA?
check <- str_subset(ratio_cols, "X")
check #6 combinations left 

# get vector of all non-ratio cols
a <- names(df)
non_ratio_cols <- a[!a %in% ratio_cols]

df_non_ratio <- df %>%
  dplyr::select(non_ratio_cols)
# get all ratio cols
df_ratio_cols <- df %>%
  dplyr::select(ratio_cols)


#remove redundent col rows - these are rows where call the remaining cols of interst are NA
#code below filters to rows where at least one of the relevent cols have data
df_filtered <- df %>% 
  dplyr::filter_at(vars(ratio_cols), any_vars(!is.na(.)))



# create dataframes with 
ratio_cols_114_115 <- str_subset(ratio_cols, pattern = "114.115")
ratio_cols_114_116 <- str_subset(ratio_cols, pattern = "114.116")
ratio_cols_114_117 <- str_subset(ratio_cols, pattern = "114.117")
ratio_cols_115_116 <- str_subset(ratio_cols, pattern = "115.116")
ratio_cols_115_117 <- str_subset(ratio_cols, pattern = "115.117")
ratio_cols_116_117 <- str_subset(ratio_cols, pattern = "116.117")
# arrange ratio_cols into a list 
ratio_cols_list <- ls(pattern = "ratio_cols_")
ratio_cols_list <- mget(ratio_cols_list)
# create a list of 6 dataframes that only conatin 1 ratio 
df_list <- map(ratio_cols_list, ~cbind(df_non_ratio, dplyr::select(df_ratio_cols, .x)))
# remove rows where the ratio columns are NA 
df_list <- map2(df_list, ratio_cols_list, ~ drop_na(.x, .y))


map(df_list, ~ names(.x))

# figure out how to adapt code above to work with list below at some point
# list(df_filtered) %>%
#   rep(length(dfs_to_list)) %>%
#   enframe(name = "id", value = "df_filtered") %>%
#   mutate(labels = dfs_to_list) -> df_filtered_list
# 
# df_filtered_list
# 
# list(df_list) %>%
#   rep(length(dfs_to_list)) %>%
#   enframe(name = "id", value = "df_list") %>%
#   mutate(labels = dfs_to_list) -> test
# test

# remove cols that only contain NAs
# map(df_list, ~ names(.x))
# 
# for (n in 1:length(dfs_to_list)) {
#   df_list[[n]] <- select_if(df_list[[n]], not_all_na)
# }
# 
# map(df_list, ~ names(.x))



# below shows Pval filter is working
a <- df_list$ratio_cols_114_115 %>%
  dplyr::filter(PVal.114.115 < 0.05)

test <- dplyr::filter(df_list$ratio_cols_114_115, df_list$ratio_cols_114_115[,13] < 0.05)

b <- df_list$df_114_vs_115 %>%
  dplyr::select(contains("PVal"))
# get list of PVal col names
pval_cols_list <- map(df_list, ~names(dplyr::select(.x, contains("PVal"))))
# note: the rlang::sym function converts strings to symbols for filter
pval_cols_list <- map(pval_cols_list, ~rlang::sym(.x))

# pval_cols_list <- c("PVal.114.115", "PVal.114.116", "PVal.114.117", "PVal.115.116", "PVal.115.117", "PVal.116.117")
# pval_cols_list <- map_chr(pval_cols_list, ~rlang::sym(.x))
# 
# substr(pval_cols_list[1], 2, 11)

# filter to rows with pval below 0.05
df_list_filtered2 <- map2(df_list, pval_cols_list, ~ dplyr::filter(.x, .y < 0.05)) # why doesn't this work?
# pval cols is 13th column, the following filters by that
df_list_filtered <- map(df_list, ~ dplyr::filter(.x, .x[,13] < 0.05))


# arrange data into tibble with a label column 
data_tibble <- tibble(df_list_filtered)
data_tibble
# check col names 
map(data_tibble$df_list_filtered, ~names(.x))
# add labels to tibble 
data_tibble$label <- c("comparing_114-115", "comparing_114-116", "comparing_114-117", 
                       "comparing_115-116", "comparing_115-117", "comparing_116-117")

# rename tibble 
run2_data_tibble <- data_tibble



# convert ratios below 1 to fold change 
# make a copy of ratio column 
run2_data_tibble$fold_change <- map(run2_data_tibble$df_list_filtered, ~abs(.x[,12]))
run2_data_tibble$fold_change[[1]]
# 1 divided by values less than 1 to get fold change 
run2_data_tibble$fold_change <- map(run2_data_tibble$fold_change, ~ ifelse(.x < 1, 1/.x, .x))
run2_data_tibble$fold_change[[1]]
# cbind copy to rest of dataframe
run2_data_tibble$df_with_fold_change <- map2(run2_data_tibble$df_list_filtered, run2_data_tibble$fold_change, ~ cbind(.x, .y))

# make indicator column for transfomation 
# create copy colum of ratio col
run2_data_tibble$df_with_fold_change <- map(run2_data_tibble$df_with_fold_change, ~ dplyr::mutate(.x, negative_fold_change = .x[,12]))
run2_data_tibble$df_list_filtered[[1]]
# create indicator of transformation
run2_data_tibble$fold_change_indicator <- map(run2_data_tibble$df_with_fold_change, ~ ifelse(.x$negative_fold_change < 1, 1, 0))
run2_data_tibble$fold_change_indicator
# remove copy col 
run2_data_tibble$df_with_fold_change <- lapply(run2_data_tibble$df_with_fold_change, 
                                               function(x) x[!(names(x) %in% c("negative_fold_change"))])
run2_data_tibble$df_with_fold_change[[1]]

# rename transformed column
run2_data_tibble$df_with_fold_change <- map(run2_data_tibble$df_with_fold_change, ~ dplyr::rename(.x, fold_change = .y))
run2_data_tibble$df_with_fold_change[[1]]

# cbind indicator column
run2_data_tibble$df_with_fold_change <- map2(run2_data_tibble$df_with_fold_change, run2_data_tibble$fold_change_indicator, 
                                             ~cbind(.x, .y))
run2_data_tibble$df_with_fold_change[[1]]
# rename indicator column 
run2_data_tibble$df_with_fold_change <- map(run2_data_tibble$df_with_fold_change, ~dplyr::rename(.x, negative_fold_change_indicator = .y))
run2_data_tibble$df_with_fold_change[[1]]

# filter out fold_changes less than 1.25, these are too small to be validated in ELISAs
run2_data_tibble$df_with_fold_change_filtered <- map(run2_data_tibble$df_with_fold_change, ~ dplyr::filter(.x, fold_change > 1.25))
run2_data_tibble

# convert negative fold changes to negative values 
run2_data_tibble$df_with_fold_change_filtered <- map(run2_data_tibble$df_with_fold_change_filtered, ~ dplyr::mutate(.x, negative_fold_change = fold_change))
negative_fold_change <- map(run2_data_tibble$df_with_fold_change_filtered, 
                            ~ ifelse(.x$negative_fold_change_indicator == 1, .x$negative_fold_change * -1, .x$negative_fold_change))
negative_fold_change
# remove negative fold change col 
run2_data_tibble$df_with_fold_change_filtered <- lapply(run2_data_tibble$df_with_fold_change_filtered, 
                                                        function(x) x[!(names(x) %in% c("negative_fold_change"))])
run2_data_tibble$df_with_fold_change_filtered[[1]]
# cbind correct negative fold change back in 
run2_data_tibble$df_with_fold_change_filtered <- map2(run2_data_tibble$df_with_fold_change_filtered, negative_fold_change, 
                                                      ~cbind(.x, .y))
run2_data_tibble$df_with_fold_change_filtered[[1]]
# rename col 
run2_data_tibble$df_with_fold_change_filtered <- map(run2_data_tibble$df_with_fold_change_filtered, 
                                                     ~dplyr::rename(.x, negative_fold_change = .y))

# remove uneeded tibble cols 
run2_data_tibble$fold_change <- NULL
run2_data_tibble$fold_change_indicator <- NULL
run2_data_tibble

# save tibble for analysis 
save(run2_data_tibble, file = "../proteomics/data/run2_data_tibble.rda")
