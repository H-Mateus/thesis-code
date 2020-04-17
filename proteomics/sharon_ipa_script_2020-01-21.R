# protein labels to network plot

# bioconductor installation - only run if it isn't installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")
# BiocManager::install("STRINGdb")

#load packages
library(tidyverse)
library(purrr)
library(STRINGdb)


# load data - NOTE: change this path to where the spreadsheet of the data is located, including the file names itself with ".xlsx"
sharon_ipa_data <- readxl::read_xlsx(path = "/run/media/mateus/14D1-5FB9/proteomic_ipa_network_plot_r/Mateus Type I, II & III data for IPA.xlsx")

# fix column names
names(sharon_ipa_data) # check before 
names(sharon_ipa_data) <- gsub(" ", "_", names(sharon_ipa_data)) # replace spaces with underscores
names(sharon_ipa_data) <- tolower(names(sharon_ipa_data)) # make all lowercase 
names(sharon_ipa_data) <- gsub("-", "_", names(sharon_ipa_data)) # replace "-" with "_"
names(sharon_ipa_data) # check after

# The function below filters rows by p-values < 0.05, converts fold changes less than
# 1 to 1/themselves, then filters fold changes less than 1.25 and finally
# converts those fold changes to negative values

# Note: to use the following function insert the name of the dataframe,
# sharon_ipa_data in this case, and the column of interest as given by the above
# names() call, type_i_fold_change and type_i_p_value respectively. 
sharon_data_filter_function <- function(dataframe = NULL, fold_change_column = NULL, p_value_column = NULL) {
  
  #print(head(dataframe$p_value_column))
  
  #p_value_column <- paste(dataframe, "$", p_value_column, sep = "")
  
  #print(head(p_value_column))
  
  # filter out type 1 p-vales less than 0.05
  sharon_ipa_data_filtered <- dplyr::filter(dataframe, p_value_column <= 0.05)
  
  #print(head(sharon_ipa_data_filtered))
  
  # Note: using square-bracket notation with a tibble results in a different output to a dataframe!
  sharon_ipa_data_filtered <- as.data.frame(sharon_ipa_data_filtered)
  
  #fold_change_column <- paste("sharon_ipa_data_filtered$", fold_change_column, sep = "")
  
  # create indicator column for rows where the type 1 fold change is less than 1
  sharon_ipa_data_filtered$fold_change_indicator <- ifelse(sharon_ipa_data_filtered[,fold_change_column] < 1, 1, 0) 
  
  #print(head(sharon_ipa_data_filtered))
  #print(str(sharon_ipa_data_filtered))
  
  # rename fold_change column to the following code works 
  sharon_ipa_data_filtered$fold_change_col <- as.numeric(sharon_ipa_data_filtered[,fold_change_column])
  
  #print(head(sharon_ipa_data_filtered))
  #print(str(sharon_ipa_data_filtered))
  #print(head(sharon_ipa_data_filtered))
  
  # convert fold change less than 1 to negative fold change
  sharon_ipa_data_filtered$fold_change <- ifelse(sharon_ipa_data_filtered$fold_change_indicator == 1, 1 / sharon_ipa_data_filtered$fold_change_col, 
                                                 sharon_ipa_data_filtered$fold_change_col)
  
  #print(head(sharon_ipa_data_filtered))
  #print(str(sharon_ipa_data_filtered))
  #print(sharon_ipa_data_filtered$fold_change)
  
  # filter out fold changes less than 1.25
  sharon_ipa_data_filtered <- sharon_ipa_data_filtered %>%
    dplyr::filter(fold_change > 1.25)
  
  # multiply negative fold changes by -1
  sharon_ipa_data_filtered$fold_change <- ifelse(sharon_ipa_data_filtered$fold_change_indicator == 1, sharon_ipa_data_filtered$fold_change * -1, 
                                                       sharon_ipa_data_filtered$fold_change)
  
  # convert to dataframe as string_db won't work otherwise
  sharon_ipa_data_filtered <- as.data.frame(sharon_ipa_data_filtered)
  
  return(sharon_ipa_data_filtered)
}


# below is how to use the function for the type I group #
names(sharon_ipa_data)
sharon_ipa_hits_type1 <- sharon_data_filter_function(dataframe = sharon_ipa_data, fold_change_column = "type_i_fold_change", 
                         p_value_column = sharon_ipa_data$type_i_p_value) 
# below is type 2
sharon_ipa_hits_type2 <- sharon_data_filter_function(dataframe = sharon_ipa_data, fold_change_column = "type_ii_fold_change", 
                                                  p_value_column = sharon_ipa_data$type_ii_p_value)
# below is type 3
sharon_ipa_hits_type3 <- sharon_data_filter_function(dataframe = sharon_ipa_data, fold_change_column = "type_iii_fold_change", 
                                                     p_value_column = sharon_ipa_data$type_iii_p_value)



# below is setting the species by NCBI taxonomy (9606 is human), and the score
# threshold - should look into what an appropriate value there is
string_db <- STRINGdb$new( version="10", species=9606,
                           score_threshold=0, input_directory="" )

# get STRING database labels - Note: this code takes a few seconds to run, the
# more proteins the longer it will take. It also gives a warning for any genes
# it coudln't match
sharon_ipa_hits_type1 <- string_db$map(sharon_ipa_hits_type1, "protein_accession")
sharon_ipa_hits_type2 <- string_db$map(sharon_ipa_hits_type2, "protein_accession")
sharon_ipa_hits_type3 <- string_db$map(sharon_ipa_hits_type3, "protein_accession")


# get network plot in r - this takes a while to run 
#string_db$plot_network(sharon_ipa_hits_type1$STRING_id)
#string_db$plot_network(sharon_ipa_hits_type2$STRING_id)
#string_db$plot_network(sharon_ipa_hits_type3$STRING_id)

# save a png of network plot - NOTE: the file path must be changed to where you
# want the .png to be saved. Don't forget to add the name of the image itself
# and add ".png" to the end
string_db$get_png(sharon_ipa_hits_type1$STRING_id, file = "/run/media/mateus/14D1-5FB9/proteomic_ipa_network_plot_r/sharon_type_i_network_plot.png")
string_db$get_png(sharon_ipa_hits_type2$STRING_id, file = "/run/media/mateus/14D1-5FB9/proteomic_ipa_network_plot_r/sharon_type_ii_network_plot.png")
string_db$get_png(sharon_ipa_hits_type3$STRING_id, file = "/run/media/mateus/14D1-5FB9/proteomic_ipa_network_plot_r/sharon_type_iii_network_plot.png")

# get labels with color of up- or down-regulated proteins
sharon_ipa_hits_type1 <- string_db$add_diff_exp_color(sharon_ipa_hits_type1, logFcColStr = "fold_change")
sharon_ipa_hits_type2 <- string_db$add_diff_exp_color(sharon_ipa_hits_type2, logFcColStr = "fold_change")
sharon_ipa_hits_type3 <- string_db$add_diff_exp_color(sharon_ipa_hits_type3, logFcColStr = "fold_change")

# save the data as .rda and .csv
save(sharon_ipa_hits_type1, file = "/run/media/mateus/14D1-5FB9/proteomic_ipa_network_plot_r/type_1_data_filtered.rda")
save(sharon_ipa_hits_type2, file = "/run/media/mateus/14D1-5FB9/proteomic_ipa_network_plot_r/type_2_data_filtered.rda")
save(sharon_ipa_hits_type3, file = "/run/media/mateus/14D1-5FB9/proteomic_ipa_network_plot_r/type_3_data_filtered.rda")

write_csv(sharon_ipa_hits_type1, path = "/run/media/mateus/14D1-5FB9/proteomic_ipa_network_plot_r/type_1_data_filtered.csv")
write_csv(sharon_ipa_hits_type2, path = "/run/media/mateus/14D1-5FB9/proteomic_ipa_network_plot_r/type_2_data_filtered.csv")
write_csv(sharon_ipa_hits_type3, path = "/run/media/mateus/14D1-5FB9/proteomic_ipa_network_plot_r/type_3_data_filtered.csv")

# post payload info to the STRING server
payload_id_1 <- string_db$post_payload(sharon_ipa_hits_type1$STRING_id, colors = sharon_ipa_hits_type1$color)
payload_id_2 <- string_db$post_payload(sharon_ipa_hits_type2$STRING_id, colors = sharon_ipa_hits_type2$color)
payload_id_3 <- string_db$post_payload(sharon_ipa_hits_type3$STRING_id, colors = sharon_ipa_hits_type3$color)

# display a STRING network png with the "halo" - NOTE: green halo means gene is down-regulated and red is up-regulated
string_db$plot_network(sharon_ipa_hits_type1$STRING_id, payload_id = payload_id_1)
string_db$plot_network(sharon_ipa_hits_type2$STRING_id, payload_id = payload_id_2)
string_db$plot_network(sharon_ipa_hits_type3$STRING_id, payload_id = payload_id_3)



# make a bar chart of fold change
ggplot(sharon_ipa_hits_type2, aes(x = protein_accession, y = fold_change)) +
  geom_bar(stat="identity", position=position_dodge()) 


# clustering of protein interaction network
cluster_list_1 <- string_db$get_clusters(sharon_ipa_hits_type1$STRING_id)
map(cluster_list_1, ~string_db$plot_network(.x))











##### Uniprot function scrap ########

# scrape part of a web page
# function to get protein description 
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
# function to get protein name 
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

# create dataframe to save results 
protein_descriptions <- data.frame("description" = c(rep("a", nrow(sharon_ipa_hits_type1))))

# get protein description - note: 
protein_descriptions$description <- map(sharon_ipa_hits_type1$protein_accession, ~ find_protein_description(.x))
# get protein name
protein_descriptions$name <- map(sharon_ipa_hits_type1$protein_accession, ~ find_protein_name(.x))
# give the list of protein functions the protein name 
protein_descriptions$description <- map2(protein_descriptions$description, protein_descriptions$name, ~ set_names(.x, .y))

