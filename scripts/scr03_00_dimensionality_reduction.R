#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 4. This script carries
# out dimensionality reduction.
# PUBLICATION: https://f1000research.com/articles/6-748

# Read command line arguments passed to the script. 
cmd_args = commandArgs(trailingOnly=TRUE)

# Temp. If using RStudio.
# setwd("/Users/farmerr2/locus/sandbox/projects/mini")


# Read yaml file. 
suppressMessages(library(yaml))

# For interactive mode
yaml_file = "gum_20.yaml"

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
tsne_no_cells <- yam$tsne_no_cells
umap_no_cells <- yam$umap_no_cells

# Load external libraries.
cat("Loading external libraries.\n")
suppressMessages(library(cytofWorkflow))
suppressMessages(library(tidyverse))

# Assigning results folders (make sure that the script1 is executed before.)
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

# Load daFrame with clustering results saved in script2.
print("Loading SCE with clustering results saved in script2.")
sce <- readRDS(file.path(results_folder, "sce_clust.rds"))

# run t-SNE & UMAP                           
set.seed(1234)          
print("Calculating t-SNE and UMAP.")
sce_tsne <- runDR(sce, "TSNE", cells = tsne_no_cells, features = "type") 
sce_umap <- runDR(sce, "UMAP", cells = umap_no_cells, features = "type") 

# Save SCE with DR data. 
print("Saving SCE with DR data.")
saveRDS(sce_tsne, file.path(results_folder, "sce_tsne.rds"))
saveRDS(sce_umap, file.path(results_folder, "sce_umap.rds"))

print("Done")
