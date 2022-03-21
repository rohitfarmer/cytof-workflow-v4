#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 4. This script carries out clustering
# using flowSOM and culterConensusPlus. 
# PUBLICATION: https://f1000research.com/articles/6-748

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

# Load daFrame with arcsinh transformed data saved in script1.
print("Loading SCE with arcsinh transformed data saved in script1.")
sce_arcsinh <- readRDS(file.path(results_folder, "sce_arcsinh.rds"))

# CLUSTERING (FlowSom).
print("CLUSTERING (FlowSOM)")
sce_clust <- cluster(sce_arcsinh, features = type_markers(sce_arcsinh), 
               xdim = 10, ydim = 10, maxK = no_of_clusters, seed = 1234) 

# Save daFrame 
print("Saving SCE with clustering data.")
saveRDS(sce_clust, file.path(results_folder, "sce_clust.rds"))

print("Done")

