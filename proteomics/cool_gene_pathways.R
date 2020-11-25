library(KEGGREST) # for pathways
library(org.Hs.eg.db) # for human ids
library(tidyverse) # for wrangling etc
library(DT) # to make nice interactive tables visible in browser
library(png) # to write png files

#########
# set working dir (use a separate folder as many images will be generated)
# make sure the entrez_id_run1_2020-02-12.rda is in the folder
#setwd("/Users/paulcool/Desktop/proteins") # this may need changing depending of file location
#########

# pathways with their gene IDs
# hsa is homo sapiens (almost 32,000 genes)
hsa_pathways  = keggLink("pathway", "hsa") %>% 
    tibble(pathway = ., id = sub("hsa:", "", names(.)))
# have a look at the pathways (above takes a little time)
hsa_pathways

# add gen Symbol and Ensemble identifiers to pathway
hsa_added = hsa_pathways %>%
    mutate(symbol = mapIds(org.Hs.eg.db, id, "SYMBOL", "ENTREZID"),
        ensembl = mapIds(org.Hs.eg.db, id, "ENSEMBL", "ENTREZID"))
# have a look:
hsa_added

# pathway descriptions (over 300 in database)
hsa_pathway_descriptions = keggList("pathway", "hsa") %>% 
    tibble(pathway = names(.), description = .)
hsa_pathway_descriptions

# now cross reference (link) pathway descriptions with the gene identifiers:
hsa_annotated = left_join(hsa_added, hsa_pathway_descriptions)
hsa_annotated

# create an interactive table in browser
datatable(hsa_annotated)

# find the unique proteins in the assays:
load("../proteomics/data/entrez_id_run1_2020-02-12.rda")
protein_ids = entrez_id_run1
protein_ids = unlist(protein_ids)
unique_proteins = unique(protein_ids)
length(unique_proteins)

# cross reference with the pathways
hsa_annotated_filtered = hsa_annotated[hsa_annotated$id %in% unique_proteins, ]

# select unique pathways
unique_pathways = hsa_annotated_filtered %>%
	dplyr::select(pathway) %>%
	distinct(pathway)
unique_pathways

# set directoy to images to be saved
setwd("../proteomics/figures")
# find the pathways (keggGet function only does max of blocks of 10)
# so map or lapply not great
# so for next loop and save images in the folder
for (n in 1:nrow(unique_pathways)){
	cat(n)
	image = keggGet(unique_pathways[n,], "image")
	writePNG(image, paste(n, "image.png", sep = ""))
	}


	
	