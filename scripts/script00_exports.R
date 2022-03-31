#!/usr/bin/env Rscript --vanilla

# PURPOSE: Export data from a single cell experiment object 
# to a data frame. This script is not a part of the Robinson's 
# pipeline.

# Read command line arguments passed to the script. 
cmd_args = commandArgs(trailingOnly=TRUE)

# Read yaml file. 
suppressMessages(library(yaml))

# For interactive mode
yaml_file = "pheno_covid_flu_all_gender.yaml"

if(interactive()){
        cat("Running in interactive mode.\n")
        yam <- read_yaml(file.path("meta", yaml_file), fileEncoding = "UTF-8") # Change yaml file for interactive execution.
}else{
        cat("Running in Rscript mode.\n")
        if (length(cmd_args) < 1){
        cat("Missing command line argument(s).\n")
        cat("Usage: script.R name.yaml\n")
        stopifnot(length(cmd_args) > 1)
        }else{
                yam <- read_yaml(cmd_args[1], fileEncoding = "UTF-8") # Pass yaml file as a command line argument.
        }
}

analysis_name <- yam$analysis_name
no_of_clusters <- yam$no_of_clusters
meta_string <- yam$meta_string

# Load external libraries.
cat("Loading external libraries.\n")
suppressMessages(library(cytofWorkflow))
suppressMessages(library(tidyverse))

# Assigning results folders (make sure that the script1 is executed before.)
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

sce_export <- function(sce, meta_string, som_string, merging_string){
        # Fetch data from SCE object.
        print("Exporting data from SCE to a tibble.")
        samp_som <- colData(sce) %>%
                as_tibble()
        clust_code <- cluster_codes(sce) %>%
                as_tibble()
        exprs_t <- t(assays(sce)$exprs) %>%
                as_tibble()
        colnames(exprs_t) <- str_replace_all(colnames(exprs_t), "-", "_")

        dat <- bind_cols(samp_som, exprs_t)
        if(is.null(merging_string)){
                som_meta <- dplyr::select(clust_code, all_of(som_string), all_of(meta_string))
        }
        else {
                som_meta <- dplyr::select(clust_code, all_of(som_string), all_of(meta_string), all_of(merging_string))
        }
        dat_meta <- left_join(dat, som_meta, by  = c("cluster_id" = som_string))

        return(dat_meta)
}

# Load post clustering data.
sce <- readRDS(file.path(results_folder, "sce_clust.rds"))

#sce_tibble <- sce_export(sce, meta_string, "som100", merging_string = "merging1")
sce_tibble <- sce_export(sce, meta_string, "som100", merging_string = NULL)

print("Saving exported data.")
saveRDS(sce_tibble, file.path(results_folder, "sce_post_clustering_export.rds"))

