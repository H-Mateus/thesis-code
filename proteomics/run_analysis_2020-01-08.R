# load packages 
library(tidyverse)
library(lubridate)
library(purrr)

# load data
load(file = "../proteomics/data/run1_data_tibble_2020-01-21.rda")
load(file = "../proteomics/data/run2_data_tibble_2020-01-21.rda")


# source custom ggplot theme
source(file = "scripts/custom_ggplot_theme_2020-01-08.R")

# # example of bar plot 
# ggplot(run1_data_tibble$df_list_filtered[[1]], aes(x = Accession, y = X114.115)) +
#   geom_bar(stat="identity", position=position_dodge())

# create bar graph of all ratios 
run1_ratio_plots <- map2(run1_data_tibble$df_with_fold_change_filtered, run1_data_tibble$label, 
                         ~ ggplot(.x ,aes(x = Accession, y = negative_fold_change)) +
                              geom_bar(stat="identity", position = position_dodge()) +
                              ggtitle(paste("Run 1", .y, sep = " ")))

run1_ratio_plots

run2_ratio_plots <- map2(run2_data_tibble$df_with_fold_change_filtered, run1_data_tibble$label, 
                         ~ ggplot(.x ,aes(x = Accession, y = negative_fold_change)) +
                           geom_bar(stat="identity", position = position_dodge()) +
                           ggtitle(paste("Run 2", .y, sep = " ")))

run2_ratio_plots


# load protein function and name dfs
load("../proteomics/data/protein_functions_run1_2020-01-20.rda")
load("../proteomics/data/protein_functions_run2_2020-01-20.rda")

# to get the first comparison
protein_functions_run1[1,1]


# test unlisting one row of data
x <- protein_functions_run1[1,]
# add as many rows as there are elements in the list 
x[nrow(x)+rep(1:(length(x$description[[1]]) - 1)),] <- NA
# unlist each column 
x[,] <- lapply(x[,], unlist)
# rep grouplabel to all rows
x$group_labels <- rep(x$group_labels[1], nrow(x))
x

# function for converting group_label column to list with labels repeated as
# many times as there are rows in the repective group
create_group_label_list <- function(dataframe) {
  
  require(purrr)
  
  group_label_list <- as.list(dataframe$group_labels)
  # get list of lengths for each group
  group_list_length <- map(dataframe$description, ~ length(.x))
  # get sum of list length 
  group_list_length_sum <- Reduce("+", group_list_length)
  # create a list with the group label reped as many times as there are elements in the respective group
  group_label_list <- map2(group_label_list, group_list_length, ~ rep(.x, .y))
  
  # add new group label list column 
  dataframe$group_labels <- group_label_list
  return(dataframe)
}

protein_functions_run1 <- create_group_label_list(protein_functions_run1)
protein_functions_run2 <- create_group_label_list(protein_functions_run2)
# rbind dataframe together
protein_descriptions <- rbind(protein_functions_run1, protein_functions_run2)

# get the sum of the list lengths
group_list_length <- map(protein_descriptions$description, ~ length(.x))
group_list_length_sum <- Reduce("+", group_list_length)
# add as many rows as there are elements in the list - end nrow() should be equal to group_list_length_sum
protein_descriptions[nrow(protein_descriptions)+rep(1:(group_list_length_sum - nrow(protein_descriptions))),] <- NA
protein_descriptions[,] <- lapply(protein_descriptions[,], unlist)

# add count of protein occurance 
protein_descriptions <- protein_descriptions %>%
  dplyr::group_by(protein_name) %>%
  dplyr::mutate(protein_count = n())

# get the dataframe in a wide format so each protein only occurs once
protein_descriptions_wide <- protein_descriptions %>%
  pivot_wider(names_from = group_labels, values_from = fold_change)

# save the long and wide formates
write_csv(protein_descriptions_wide, path = "../proteomics/data/protein_descriptions_wide.csv")
save(protein_descriptions, file = "../proteomics/data/protein_descriptions_long.rda")
save(protein_descriptions_wide, file = "../proteomics/data/protein_descriptions_wide.rda")

# make a plot of fold changes
protein_descriptions %>%
  dplyr::group_by(group_labels) %>%
  ggplot(aes(x = protein_name, y = fold_change)) + 
  geom_bar(stat="identity", position = position_dodge())

protein_descriptions$group_labels <- as.factor(protein_descriptions$group_labels)

# make data ordered by increasing fold change 
#protein_descriptions$protein_name <- factor(protein_descriptions$protein_name, 
 #                           levels = protein_descriptions$protein_name[order(protein_descriptions$fold_change)])

# make list of bar graphs for each group
plot_list <- lapply(sort(unique(protein_descriptions$group_labels)), function(i) {
  ggplot(protein_descriptions[protein_descriptions$group_labels == i,], aes(x = reorder(protein_name, - fold_change), y = fold_change)) + 
    geom_bar(stat="identity", position = position_dodge())

})

plot_list



#data3$x <- factor(data3$x,                                    # Factor levels in decreasing order
 #                 levels = data3$x[order(data3$y, decreasing = TRUE)])

#ggplot(data3, aes(x, y)) +                                    # Decreasi




# STRINGdb code 

# map gene names to STRING database 
library(STRINGdb)

# below is setting the species by NCBI taxonomy (9606 is human), and the score
# threshold - should look into what an appropriate value there is
string_db <- STRINGdb$new( version="10", species=9606,
                            score_threshold=0, input_directory="" )

# example use of STRINGdb below
# map gene labels to STRING labels
test_map <- string_db$map(run1_data_tibble$df_with_fold_change_filtered[[1]], "gene_name")
# extract STRING labels
test_hits <- test_map$STRING_id
# create network plot
test_plot <- string_db$plot_network(test_hits)

#test <- string_db$get_graph()
string_db$get_png(test_hits, file = "test_png_title.png")

ll <- list.files(patt='*.png')


add_file_name_to_png <- function(file_name = NULL, file_directory = NULL){
  require(png)
  require(grid)
  setwd(file_directory)
  img <- as.raster(readPNG(file_name))
  ## get the file name
  x_name <- gsub('(.*).png','\\1',file_name)
  ## new device for new image version
  png(file =paste(x_name, '_modified', '.png', sep=''))
  grid.raster(img)
  ## here I add title - cex determins text size
  grid.text(label = x_name, x = 0.5, y = 0.9, gp = gpar(cex = 2))
  dev.off()
  
}

add_file_name_to_png(file_name = "test_png_title.png")

file_vector <- c("test_png_title.png", "test_png_title2.png")

imgs <- map(file_vector, ~add_file_name_to_png(.x))


string_db$get_link(test_hits) ### NOTE: this is the function to get the url


# get labels with color of up- or down-regulated proteins
example1_mapped_pval <- string_db$add_diff_exp_color(test_map, logFcColStr = "negative_fold_change")
# post payload info to the STRING server
payload_id <- string_db$post_payload(example1_mapped_pval$STRING_id, colors = example1_mapped_pval$color)
# display a STRING network png with the "halo" - NOTE: green halo means gene is down-regulated and red is up-regulated
string_db$plot_network(test_hits, payload_id = payload_id)
# clustering
cluster_list <- string_db$get_clusters(run1_data_tibble$protein_interactions[[1]]$STRING_id)
map(cluster_list, ~string_db$plot_network(.x))



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

# get vector of file names
figure_file_vector <- unlist(list.files(path = "../proteomics/figures/", pattern = "*.png"))
# add file name to pngs
imgs <- map(figure_file_vector, ~ add_file_name_to_png(file_name = .x, file_directory = "/home/mateus/git/mateus/proteomics/figures/"))
# return working directory 
setwd("../../500_patient/")





# make network plots
run1_network_plots <- map(run1_hits, ~ string_db$plot_network(.x))

# add up-/down-regulation colors to network plot 
run1_data_tibble$protein_interactions_color <- map(run1_data_tibble$protein_interactions, ~ 
                                                     string_db$add_diff_exp_color(.x, logFcColStr = "negative_fold_change"))
run2_data_tibble$protein_interactions_color <- map(run2_data_tibble$protein_interactions, ~ 
                                                     string_db$add_diff_exp_color(.x, logFcColStr = "negative_fold_change"))
# post payload inforamtion to the STRING server
run1_data_tibble$payload_id <- map(run1_data_tibble$protein_interactions_color, ~ 
                                     string_db$post_payload(.x$STRING_id, colors = .x$color))
run2_data_tibble$payload_id <- map(run2_data_tibble$protein_interactions_color, ~ 
                                     string_db$post_payload(.x$STRING_id, colors = .x$color))

# get list of urls for netowrk plots - NOTE: get_link only works correctly with
# a character vector and run1_hits is a list of dataframes - also note that get link doesn't work with NAs
run1_hits_no_na <- map(run1_hits, ~ .x[!is.na(.x)])
run2_hits_no_na <- map(run2_hits, ~ .x[!is.na(.x)])

network_urls_run1 <- map2(run1_hits_no_na, run1_data_tibble$payload_id, ~ string_db$get_link(.x, payload_id = .y))
network_urls_run2 <- map2(run2_hits_no_na, run2_data_tibble$payload_id, ~ string_db$get_link(.x, payload_id = .y))



### test to download scallable image from STRING ###
library(rvest)
# get url 
url <- network_urls_run1[[2]]

# function to take a url from string_db and save an svg of the protein
# interaction network plot
get_stringdb_svg <- function(url, file_path) {
  # get packages 
  require(rvest)
  require(data.table)
  # read html from the website 
  webpage <- read_html(url)
  # select svg css element
  svg_css <- html_nodes(webpage, ".row:nth-child(3) a") 
  # select the attribute from the node - this gives the download link
  svg_link <- html_attr(svg_css, "href")
  # add begining of web address
  svg_full_link <- paste("version10.string-db.org", svg_link, sep = "")
  # downlead the svg from the link 
  download.file(svg_full_link, file_path)
}

# make list of file paths
svg_file_paths_run1 <- as.list(run1_data_tibble$label)

# use svg function 
map2(network_urls_run1, file_name_list_run1, ~ get_stringdb_svg(.x, paste("../proteomics/figures/network_svgs/", .y, today(), ".svg", sep = "_")))
map2(network_urls_run2, file_name_list_run2, ~ get_stringdb_svg(.x, paste("../proteomics/figures/network_svgs/", .y, today(), ".svg", sep = "_")))


# display network plot of "halo"
run1_network_plots_colored <- map2(run1_hits, run1_data_tibble$payload_id, ~string_db$plot_network(.x, payload_id = .y))

# cluster analysis 
run1_data_tibble$cluster_list <- map(run1_data_tibble$protein_interactions, ~ string_db$get_clusters(.x$STRING_id))
run1_cluster_plots <- pmap(run1_data_tibble$cluster_list, ~ string_db$plot_network(.x))




# Reactome code for pathway analysis 

library(ReactomePA)
library(org.Hs.eg.db) # this package allows you to get entrez gene id from other label types
# see the types of labels 
columns(org.Hs.eg.db)
# get the entrez gene id from the uniprot label
c <- mapIds(org.Hs.eg.db, run1_data_tibble$protein_interactions_color[[1]]$protein_label, 'ENTREZID', 'UNIPROT')

# enrichment pathway analysis
x <- enrichPathway(gene = c)
# visualize enrichment result
barplot(x, showCategory = 8)
dotplot(x, showCategory = 15)
emapplot(x)
cnetplot(x, categorySize = 'pvalue', foldChange = run1_data_tibble$protein_interactions_color[[1]]$negative_fold_change)
# example with cnetplot
library(DOSE)
data("geneList") # note that geneList is a named numeric where the numbers are fold changes and the names are gene labels
de <- names(geneList)[1:100]
x <- enrichDO(de)
cnetplot(x, showCategory = 15, foldChange = geneList)

# do the same as above for the whole run 
# get entrez ids
entrez_id_run1 <- map(run1_data_tibble$protein_interactions_color, ~ mapIds(org.Hs.eg.db, .x$protein_label, 'ENTREZID', 'UNIPROT')) 
# get enrichment pathway
reactome_df_run1 <- map(entrez_id_run1, ~ enrichPathway(gene = .x, pAdjustMethod = "bonferroni"))
# visualize 
barplots_run1 <- map(reactome_df_run1, ~ barplot(.x, showCategory = 8))
dotplots_run1 <- map(reactome_df_run1, ~ dotplot(.x, showCategory = 15))
emapplot_run1 <- map(reactome_df_run1, ~ emapplot(.x))
cnetplot_run1 <- map2(reactome_df_run1, run1_data_tibble$protein_interactions_color, ~ cnetplot(.x, foldChange = .y$negative_fold_change))

# save data
save(entrez_id_run1, file = "../proteomics/data/entrez_id_run1_2020-02-12.rda")
save(reactome_df_run1, file = "../proteomics/data/reactome_df_run1_2020-02-12.rda")
save(run1_data_tibble, file = "../proteomics/data/run1_data_tibble_2020-02-12.rda")

gene_id_test <- c("718", "2", "721")
gene_pathway <- enrichPathway(gene_id_test, pAdjustMethod = "bonferroni")
emapplot(gene_pathway)

m = matrix(c(2, 1, 23, 10631), byrow = TRUE, ncol = 2)
chisq.test(m)
fisher.test(m)
# get entrez ids
entrez_id_run2 <- map(run2_data_tibble$protein_interactions_color, ~ mapIds(org.Hs.eg.db, .x$protein_label, 'ENTREZID', 'UNIPROT')) 
# get enrichment pathway
reactome_df_run2 <- map(entrez_id_run2, ~ enrichPathway(gene = .x))
# visualize 
barplots_run2 <- map(reactome_df_run2, ~ barplot(.x, showCategory = 8))
dotplots_run2 <- map(reactome_df_run2, ~ dotplot(.x, showCategory = 15))
emapplot_run2 <- map(reactome_df_run2, ~ emapplot(.x))
cnetplot_run2 <- map2(reactome_df_run2, run2_data_tibble$protein_interactions_color, ~ cnetplot(.x, foldChange = .y$negative_fold_change))

# set plot list names to describe group
barplots_run1 <- set_names(barplots_run1, file_name_list_run1)
dotplots_run1 <- set_names(dotplots_run1, file_name_list_run1)
emapplot_run1 <- set_names(emapplot_run1, file_name_list_run1)
barplots_run2 <- set_names(barplots_run2, file_name_list_run2)
dotplots_run2 <- set_names(dotplots_run2, file_name_list_run2)
emapplot_run2 <- set_names(emapplot_run2, file_name_list_run2)

# save plots 
save(barplots_run1, file = "../proteomics/figures/reactome/barplots_run1.rda")
save(dotplots_run1, file = "../proteomics/figures/reactome/dotplots_run1.rda")
save(emapplot_run1, file = "../proteomics/figures/reactome/emapplot_run1.rda")
save(barplots_run2, file = "../proteomics/figures/reactome/barplots_run2.rda")
save(dotplots_run2, file = "../proteomics/figures/reactome/dotplots_run2.rda")
save(emapplot_run2, file = "../proteomics/figures/reactome/emapplot_run2.rda")

# view plots 
barplots_run1
dotplots_run1
emapplot_run1

barplots_run2
dotplots_run2
emapplot_run2






# plot fold change 
# run1_fold_change_plots <- map2(run1_data_tibble$df_with_fold_change, run1_data_tibble$label, ~ ggplot(.x ,aes(x = Accession, y = .x[,17])) +
#                            geom_bar(stat="identity", position=position_dodge()) +
#                            ylab(paste( "Log10 -", names(.x[12]), sep = " ")) +
#                            scale_y_continuous(trans='log10', labels = scales::comma) +
#                            ggtitle(paste("Run 1", .y, sep = " ")))
# 
# run1_fold_change_plots
# 
# run2_fold_change_plots <- map2(run2_data_tibble$df_with_fold_change, run2_data_tibble$label, ~ ggplot(.x ,aes(x = Accession, y = .x[,17])) +
#                                  geom_bar(stat="identity", position=position_dodge()) +
#                                  ylab(paste( "Log10 -", names(.x[12]), sep = " ")) +
#                                  scale_y_continuous(trans='log10', labels = scales::comma) +
#                                  ggtitle(paste("Run 2 - fold change", .y, sep = " ")))
# 
# run2_fold_change_plots


# volcano plot?
# with(run1_data_tibble$df_list_filtered[[1]]$, plot(, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))


# bioconductor installation 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")

# check bioconductor version with:
BiocManager::version()
# check packages are up-to-date with:
BiocManager::valid()

#BiocManager::install("BioNet")
#BiocManager::install("DLBCL")

# library(BioNet)
# library(DLBCL)
# data(dataLym)
# data(interactome)
# 
# pvals <- cbind(t = dataLym$t.pval, s = dataLym$s.pval)
# rownames(pvals) <- dataLym$label
# pval <- aggrPvals(pvals, order = 2, plot = TRUE)
# 
# subnet <- subNetwork(dataLym$label, interactome)
# subnet <- rmSelfLoops(subnet)
# subnet
# 
# fb <- fitBumModel(pval, plot = TRUE)
# scores <- scoreNodes(subnet, fb, fdr = 0.001)
# 
# module <- runFastHeinz(subnet, scores)
# logFC <- dataLym$diff
# names(logFC) <- dataLym$label
# plotModule(module, scores = scores, diff.expr = logFC)
# 
# test_label <- c("APOB(338)", "SERPINA3(12)", "HP(3240)")
# 
# subnet2 <- subNetwork(test_label, interactome)
# subnet2<- rmSelfLoops(subnet2)
# subnet2


# r-bloggers code for scrapping NCBI 

# get.ppiNCBI <- function(g.n) {
#   require(XML)
#   ppi <- data.frame()
#   for(i in 1:length(g.n)){
#     o <- htmlParse(paste("http://www.ncbi.nlm.nih.gov/gene/", g.n[i], sep=''))
#     # check if interaction table exists
#     exist <- length(getNodeSet(o, "//table//th[@id='inter-prod']"))>0
#     if(exist){
#       p <- getNodeSet(o, "//table")
#       ## need to know which table is the good one
#       for(j in 1:length(p)){
#         int <- readHTMLTable(p[[j]])
#         if(colnames(int)[2]=="Interactant"){break}
#       }
#       ppi <- rbind(ppi, data.frame(egID=g.n[i], intSymbol=int$`Other Gene`))
#     }
#     # play nice! and avoid being kicked out from NCBI servers
#     Sys.sleep(1)
#   }
#   if(dim(ppi)[1]>0){
#     ppi <- unique(ppi)
#     print(paste(dim(ppi)[1], "interactions found"))
#     return(ppi)
#   } else{
#     print("No interaction found")
#   }
# }




