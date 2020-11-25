# scrape part of a web page

find_protein_description <- function(protein_id = ""){
	require(rvest)  # hadley wickham's package that can harvest html files
	website <- paste("https://www.uniprot.org/uniprot/", protein_id, sep = "") # website
	webpage <- read_html(website) # webpage read
	results <- webpage %>% html_nodes("meta") # get relevant data
	result <- as.character(results[[4]]) # 4th node is what you want
	result <- str_remove(result, "<meta content=\"") # strip beginning
	result <- str_remove(result, "\" name=\"description\">\n") # strip end
	return(result) # return result
}

find_protein_name <- function(protein_id = ""){
	require(rvest)  # hadley wickham's package that can harvest html files
	website <- paste("https://www.uniprot.org/uniprot/", protein_id, sep = "") # website
	webpage <- read_html(website) # webpage read
	name <- webpage %>% html_nodes("head")  %>% html_node("title") # get protein name
	name <- as.character(name[[1]])
	name <- str_remove(name, "<title>") # strip beginning
	name <- str_remove(name, "gene &amp; protein</title>\n") # strip end
	return(name) # return name
}

# make version of functions for use in list of dataframes
find_protein_description_list <- function(protein_label_column) {
  map(protein_label_column, ~ find_protein_description(.x))
}

find_protein_name_list <- function(protein_label_column) {
  map(protein_label_column, ~ find_protein_name(.x))
}

# load packages 
library(tidyverse)

# source data cleaning function
source(file = "../proteomics/scripts/proteinpilot_clean_function_2020-01-10.R")

# get data
run1_data_tibble <- proteinpilot_export_clean("run1")
run2_data_tibble <- proteinpilot_export_clean("run2")

# get protein description in new df
protein_functions_run1 <- data.frame("description" = c(rep("a", nrow(run1_data_tibble))))
protein_functions_run1$description <- map(run1_data_tibble$protein_interactions, ~ find_protein_description_list(.x$protein_label))
# get protein name 
protein_functions_run1$protein_name <- map(run1_data_tibble$protein_interactions, ~ find_protein_name_list(.x$protein_label))

# do the same for run2
protein_functions_run2 <- data.frame("description" = c(rep("a", nrow(run1_data_tibble))))
protein_functions_run2$description <- map(run2_data_tibble$protein_interactions, ~ find_protein_description_list(.x$protein_label))
# get protein name 
protein_functions_run2$protein_name <- map(run2_data_tibble$protein_interactions, ~ find_protein_name_list(.x$protein_label))

# create list of groups 
group_label_run1 <- c("c_improv_acute_vs_c_improv_subacute", "c_improv_acute_vs_c_nonimprov_acute", 
                 "c_improv_acute_vs_c_nonimprov_subacute", "c_improv_subacute_vs_c_nonimprov_acute",
                 "c_improv_subacute_vs_c_nonimprove_subacute", "c_nonimprov_acute_vs_c_nonimprov_subacute")

group_label_run2 <- c("c_improv_vs_c_nonimprov", "c_improv_vs_a", "c_improv_vs_d", "c_nonimprov_vs_a", 
                      "c_nonimprov_vs_d", "a_vs_d")
# set rownames to group comparisons
#rownames(protein_functions_run1) <- group_label_run1
#rownames(protein_functions_run2) <- group_label_run2

# create column of group labels 
protein_functions_run1$group_labels <- group_label_run1
protein_functions_run2$group_labels <- group_label_run2

# add fold change to dataframe
protein_functions_run1$fold_change <- map(run1_data_tibble$protein_interactions, ~ as.list(.x$negative_fold_change))
protein_functions_run2$fold_change <- map(run2_data_tibble$protein_interactions, ~ as.list(.x$negative_fold_change))

# give the list of protein descriptions the protein names
protein_functions_run1$description <- map2(protein_functions_run1$description, protein_functions_run1$protein_name, ~ set_names(.x, .y))
protein_functions_run2$description <- map2(protein_functions_run2$description, protein_functions_run2$protein_name, ~ set_names(.x, .y))

# save protein dataframes
save(protein_functions_run1, file = "../proteomics/data/protein_functions_run1_2020-01-20.rda")
save(protein_functions_run2, file = "../proteomics/data/protein_functions_run2_2020-01-20.rda")

# # check it works with tp53
# my_protein <- "P04637"
# find_protein_description(my_protein)
# find_protein_name(my_protein)
# 
# # now with a data frame
# df <- data.frame("protein" = c("P04637", "Q00987", "O75874", "P48735", "P35222", "Q01844"))
# df$description <- map(df$protein, ~find_protein_description(.x))
# df$name <- map(df$protein, ~find_protein_name(.x))
# df









