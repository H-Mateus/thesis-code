#load packages
library(tidyverse)
library(purrr)
library(STRINGdb)
library(lubridate)

#run = "run1"

# function for cleaning proteomic data - Note: the empty columns in the export between has to be removed for this to work
proteinpilot_export_clean <- function(run) {
  

  
  #load data - note: files are tab seperated and row.names must be set to null for
  #protein_summary files
  
  #proteins_run1 <- read.table(file = '../Proteomics/data/proteomics_run1_proteinsummary_proteinpilot_2019-10-25.txt', 
  #                            sep="\t", header=TRUE, row.names = NULL)
  
  # note: data from proteinpilot software contains empty cols and is a tab
  # seperated .txt file - if those are loaded in as they are the data is not
  # formated corretely - fixed by openeing in spreadsheet software, deleting empty
  # cols and saving as a .csv
  data_path <- "../proteomics/data/protein_pilot_export/"   # path to the data
  files <- dir(data_path, pattern = paste("proteomics_", as.character(run), sep = "")) # get file names
  files
  # create a data frame holding the file names
  data_list <- tibble(filename = files) %>% 
    mutate(file_contents = map(filename,                            # read files into
                               ~ read_tsv(file.path(data_path, .), col_names = TRUE)) # a new data column
    )  
  data_list
  
  
  # remove cols that only contain NAs
  #function below returns TRUE if there are only NAs in an object
  not_all_na <- function(x) any(!is.na(x))
  
  
  num_df <- length(data_list$filename)
  data_list$df_filtered <- list(rep(1:nrow(data_list)))
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
  #a <- df_list$ratio_cols_114_115 %>%
   # dplyr::filter(PVal.114.115 < 0.05)
  
  
  test <- dplyr::filter(df_list$ratio_cols_114_115, df_list$ratio_cols_114_115[,13] < 0.05)
  
  #b <- df_list$df_114_vs_115 %>%
   # dplyr::select(contains("PVal"))
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
  
  df_list_filtered
  
  #df_list_filtered3 <- map(df_list, ~ dplyr::filter(.x, .x[,grep("PVal" %in% names(.x))] < 0.05))
  
  
  # arrange data into tibble with a label column 
  data_tibble <- tibble(df_list_filtered)
  data_tibble
  # check col names 
  map(data_tibble$df_list_filtered, ~names(.x))
  # add labels to tibble 
  data_tibble$label <- c("comparing_114-115", "comparing_114-116", "comparing_114-117", 
                         "comparing_115-116", "comparing_115-117", "comparing_116-117")
  
  # rename tibble 
  run1_data_tibble <- data_tibble
  
  
  # convert ratios below 1 to fold change 
  # make a copy of ratio column 
  run1_data_tibble$fold_change <- map(run1_data_tibble$df_list_filtered, ~abs(.x[,12]))
  run1_data_tibble$fold_change[[1]]
  # 1 divided by values less than 1 to get fold change 
  run1_data_tibble$fold_change <- map(run1_data_tibble$fold_change, ~ ifelse(.x < 1, 1/.x, .x))
  run1_data_tibble$fold_change[[1]]
  # cbind copy to rest of dataframe
  run1_data_tibble$df_with_fold_change <- map2(run1_data_tibble$df_list_filtered, run1_data_tibble$fold_change, ~ cbind(.x, .y))
  
  # make indicator column for transfomation 
  # create copy colum of ratio col
  run1_data_tibble$df_with_fold_change <- map(run1_data_tibble$df_with_fold_change, ~ dplyr::mutate(.x, negative_fold_change = .x[,12]))
  run1_data_tibble$df_list_filtered[[1]]
  # create indicator of transformation
  run1_data_tibble$fold_change_indicator <- map(run1_data_tibble$df_with_fold_change, ~ ifelse(.x$negative_fold_change < 1, 1, 0))
  run1_data_tibble$fold_change_indicator
  # remove copy col 
  run1_data_tibble$df_with_fold_change <- lapply(run1_data_tibble$df_with_fold_change, 
                                                 function(x) x[!(names(x) %in% c("negative_fold_change"))])
  run1_data_tibble$df_with_fold_change[[1]]
  
  # rename transformed column
  run1_data_tibble$df_with_fold_change <- map(run1_data_tibble$df_with_fold_change, ~ dplyr::rename(.x, fold_change = .y))
  run1_data_tibble$df_with_fold_change[[1]]
  
  # cbind indicator column
  run1_data_tibble$df_with_fold_change <- map2(run1_data_tibble$df_with_fold_change, run1_data_tibble$fold_change_indicator, 
                                               ~cbind(.x, .y))
  run1_data_tibble$df_with_fold_change[[1]]
  # rename indicator column 
  run1_data_tibble$df_with_fold_change <- map(run1_data_tibble$df_with_fold_change, ~dplyr::rename(.x, negative_fold_change_indicator = .y))
  run1_data_tibble$df_with_fold_change[[1]]
  
  # filter out fold_changes less than 1.25, these are too small to be validated in ELISAs
  run1_data_tibble$df_with_fold_change_filtered <- map(run1_data_tibble$df_with_fold_change, ~ dplyr::filter(.x, fold_change > 1.25))
  run1_data_tibble
  
  # convert negative fold changes to negative values 
  run1_data_tibble$df_with_fold_change_filtered <- map(run1_data_tibble$df_with_fold_change_filtered, ~ dplyr::mutate(.x, negative_fold_change = fold_change))
  negative_fold_change <- map(run1_data_tibble$df_with_fold_change_filtered, 
                              ~ ifelse(.x$negative_fold_change_indicator == 1, .x$negative_fold_change * -1, .x$negative_fold_change))
  negative_fold_change
  # remove negative fold change col 
  run1_data_tibble$df_with_fold_change_filtered <- lapply(run1_data_tibble$df_with_fold_change_filtered, 
                                                 function(x) x[!(names(x) %in% c("negative_fold_change"))])
  run1_data_tibble$df_with_fold_change_filtered[[1]]
  # cbind correct negative fold change back in 
  run1_data_tibble$df_with_fold_change_filtered <- map2(run1_data_tibble$df_with_fold_change_filtered, negative_fold_change, 
                                                        ~cbind(.x, .y))
  run1_data_tibble$df_with_fold_change_filtered[[1]]
  # rename col 
  run1_data_tibble$df_with_fold_change_filtered <- map(run1_data_tibble$df_with_fold_change_filtered, 
                                                       ~dplyr::rename(.x, negative_fold_change = .y))
  
  # remove uneeded tibble cols 
  run1_data_tibble$fold_change <- NULL
  run1_data_tibble$fold_change_indicator <- NULL
  
  
  # code to extract string between two characters 
  try <- run1_data_tibble$df_with_fold_change_filtered[[1]]$Accession
  try
  library(qdapRegex)
  try <- rm_between(try, "|", "|", extract = TRUE)
  unlist(try)
  
  
  run1_data_tibble$df_with_fold_change_filtered <- map(run1_data_tibble$df_with_fold_change_filtered, 
                                                       ~dplyr::mutate(.x, protein_label = 
                                                                        unlist(rm_between(.x$Accession, "|", "|", extract = TRUE))))
  run1_data_tibble$df_with_fold_change_filtered[[1]]
  
  # code to extract gene 
  gene_try <- as.character(run1_data_tibble$df_with_fold_change_filtered[[1]]$Name)
  rm_between(gene_try, "GN=", "PE",  extract = TRUE)
  
  
  run1_data_tibble$df_with_fold_change_filtered <- map(run1_data_tibble$df_with_fold_change_filtered, 
                                                       ~dplyr::mutate(.x, gene_name = 
                                                                        unlist(rm_between(.x$Name, "GN=", "PE", 
                                                                                          extract = TRUE))))
  run1_data_tibble$df_with_fold_change_filtered[[1]]
  str(run1_data_tibble$df_with_fold_change_filtered[[1]])
  
  # filter out keratin
  
  # convert Name col to character - Note: the rm_between works with factors
  run1_data_tibble$df_with_fold_change_filtered <- map(run1_data_tibble$df_with_fold_change_filtered, 
                                                      ~dplyr::mutate(.x, Name = as.character(.x$Name)))
  # filter out keratin rows
  run1_data_tibble$df_with_fold_change_filtered <- map(run1_data_tibble$df_with_fold_change_filtered,
                                                       ~dplyr::filter(.x, !grepl("Keratin", Name)))
  
  
  return(run1_data_tibble)
  # save tibble for analysis 
  #save(run1_data_tibble, file = "../proteomics/data/run1_data_tibble.rda")
  
}

# get data
run1_data_tibble <- proteinpilot_export_clean("run1")
run2_data_tibble <- proteinpilot_export_clean("run2")

# set STRINGdb paramaters -below is setting the species by NCBI taxonomy (9606
# is human), and the score threshold - should look into what an appropriate
# value there is
string_db <- STRINGdb$new( version="10", species=9606,
                           score_threshold=0, input_directory="" )

# map gene labels to STRING database labels - note there were some failed matches
run1_data_tibble$protein_interactions <- map(run1_data_tibble$df_with_fold_change_filtered, ~ string_db$map(.x, "gene_name"))
run2_data_tibble$protein_interactions <- map(run2_data_tibble$df_with_fold_change_filtered, ~ string_db$map(.x, "gene_name"))


### there are duplicate genes for some reason - remove them 
nrow(run1_data_tibble$protein_interactions[[1]])
run1_data_tibble$protein_interactions <- map(run1_data_tibble$protein_interactions, ~ dplyr::distinct(.x, gene_name, .keep_all = TRUE))
run2_data_tibble$protein_interactions <- map(run2_data_tibble$protein_interactions, ~ dplyr::distinct(.x, gene_name, .keep_all = TRUE))
nrow(run1_data_tibble$protein_interactions[[1]])


# save tibbles as the code to get labels from the database takes time
save(run1_data_tibble, file = "../proteomics/data/run1_data_tibble_2020-01-21.rda")
save(run2_data_tibble, file = "../proteomics/data/run2_data_tibble_2020-01-21.rda")


# extract the STRING labels
run1_hits <- map(run1_data_tibble$protein_interactions, ~dplyr::select(.x, "STRING_id"))
run2_hits <- map(run2_data_tibble$protein_interactions, ~dplyr::select(.x, "STRING_id"))

# make a list of file names to save and later add to network .png images
file_name_list_run1 <- as.list(c("c_improv_acute_vs_c_improv_subacute", 
                                 "c_improv_acute_vs_c_nonimprov_acute", 
                                 "c_improv_acute_vs_c_nonimprov_subacute",
                                 "c_improv_subacute_vs_c_nonimprov_acute",
                                 "c_improv_subacute_vs_c_nonimprove_subacute",
                                 "c_nonimprov_acute_vs_c_nonimprov_subacute"))

file_name_list_run2 <- as.list(c("c_improv_vs_c_nonimprov", 
                                 "c_improv_vs_a", 
                                 "c_improv_vs_d",
                                 "c_nonimprov_vs_a",
                                 "c_nonimprov_vs_d",
                                 "a_vs_d"))
# save pngs of interaction networks
map2(run1_hits, file_name_list_run1, ~ string_db$get_png(.x, file = paste("../proteomics/figures/", .y, "_", today(), ".png", sep = "")))
map2(run2_hits, file_name_list_run2, ~ string_db$get_png(.x, file = paste("../proteomics/figures/", .y, "_", today(), ".png", sep = "")))

# run1_data_tibble <- proteinpilot_export_clean("run1")
# run2_data_tibble <- proteinpilot_export_clean("run2")
# save tibble for analysis 
#save(run1_data_tibble, file = "../proteomics/data/run1_data_tibble.rda")
#save(run2_data_tibble, file = "../proteomics/data/run2_data_tibble.rda")



  # # below is a different attempt at getting data into a list
  # # create a label column for rows which contain data on interest
  # cols_to_select <- str_remove(check, "X")
  # cols_to_select
  # 
  # df_filtered$label <- ""
  # 
  # for (i in 1:nrow(df_filtered)) {
  #   
  #   if (!is.na(df_filtered$X114.115[i])) {
  #     df_filtered$label[i] <- "compare_114-115"
  #   } else if (!is.na(df_filtered$X114.116[i])) {
  #      df_filtered$label[i] <- "compare_114-116"
  #   } else if (!is.na(df_filtered$X114.117[i])) {
  #      df_filtered$label[i] <- "compare_114-117"
  #   } else if (!is.na(df_filtered$X115.116[i])) {
  #      df_filtered$label[i] <- "compare_115-116"
  #   } else if (!is.na(df_filtered$X115.117[i])) {
  #      df_filtered$label[i] <- "compare_115-117"
  #   } else if (!is.na(df_filtered$X116.117[i])) {
  #      df_filtered$label[i] <- "compare_116-117"
  #   }
  #   
  # }
  # 
  # 
  # #revert to tibble with data lists for p-value filtering
  # df_list2 <- df_filtered %>%
  #   nest(data = c(-filename))
  # 
  # df_list2
  # a <- df_list2$data[[1]]
  # b <- df_list2$data[[2]]
  # c <- df_list2$data[[3]]
  # 
  # #add empty rows to tibble 
  # df_list2[nrow(df_list2)+3,] <- NA
  # df_list2
  # # add correct labels for dataframes
  # df_list2$label <- c("comparing_114-115", "comparing_115-116", "comparing_114-117", 
  #                 "comparing_114-116", "comparing_115-117", "comparing_116-117")
  # df_list2
  # 
  # 
  # df <- tibble(x = c(1, 1, 1, 2, 2, 3), y = 1:6, z = 6:1)
  # # Note that we get one row of output for each unique combination of
  # # non-nested variables
  # df %>% nest(data = c(y, z))
  # df %>% nest (data = -x)
  