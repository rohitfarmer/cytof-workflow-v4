#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow version 4. This script produces
# arcsinh transformed and scaled data for further analysis.
# PUBLICATION: https://f1000research.com/articles/6-748

# Read command line arguments passed to the script. 
cmd_args = commandArgs(trailingOnly=TRUE)

# For interactive mode
yaml_file = "pheno_covid_flu_all_gender.yaml"

# Read yaml file. 
suppressMessages(library(yaml))
cat("Running script01_archsinh_transform.R\n")
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
data_location <- yam$data_location
args_file <- yam$args_file
panel_file <- yam$panel_file
condition_levels <- yam$condition_levels
no_of_clusters <- yam$no_of_clusters
tsne_no_cells <- yam$tsne_no_cells
umap_no_cells <- yam$umap_no_cells
meta_string <- yam$meta_string
cofactor <- yam$cofactor

# Load external libraries.
cat("Loading external libraries.\n")
suppressMessages(library(cytofWorkflow))
suppressMessages(library(tidyverse))

# CREATE FOLDERS
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

if(!dir.exists(results_folder)){
  cat(sprintf("Creating %s folder", results_folder))
  dir.create(results_folder, recursive = TRUE)
} else {
  cat(sprintf("%s folder already exists. Output will be over written.", results_folder))
}

if(!dir.exists(figures_folder)){
  cat(sprintf("Creating %s folder", figures_folder))
  dir.create(figures_folder, recursive = TRUE)
} else {
  cat(sprintf("%s folder already exists. Output will be over written.", figures_folder))
}

# LOAD DATA
# Read experiment metadata.
sprintf("Loading experiment metadata: %s", args_file)
md<-suppressMessages(read_tsv(file.path("meta",args_file)))

# Specify levels for conditions & sample IDs to assure desired ordering.
md$condition <- factor(md$condition, levels = condition_levels)
md$sample_id  <- factor(md$sample_id, levels = md$sample_id[order(md$condition)]) # this is what paper suggests.
#md$file_name <- file.path("data", data_location, md$file_name)
cols_gen <- c("file_name", "sample_id", "condition")
diff_cols <-  setdiff(colnames(md), cols_gen)

# Read panel data.
sprintf("Loading panel information: %s", panel_file)
panel <- suppressMessages(read_tsv(file.path("meta", panel_file)))

# Load .fcs file into a flowset.
sprintf("Loading fcs files mentioned in the experiment metadata into a flowSet: %s", data_location)
fcs_raw <- read.flowSet(file = md$file_name, path = file.path("data", data_location),
                        transformation = FALSE, truncate_max_range = FALSE)

# Check for file names. They should match to what is in the md$file_name.
ids <- c(keyword(fcs_raw, "FILENAME"))
sprintf("Checking .fcs filenames in flowSet.")
sprintf(ids)

# Spot check that all panel columns are in the flowSet object
sprintf("Spot check that all panel columns are in the flowSet object %s", 
        all(panel$fcs_colname %in% colnames(fcs_raw)))

# Construct a SingleCellExperiment object. 
sce <- prepData(fcs_raw, panel, md, features = panel$fcs_colname, transform = TRUE, cofactor = cofactor,
                md_cols = list(file = "file_name", id = "sample_id", 
                               factors = c("condition", diff_cols)))

# Remove fcs raw to save memory.
rm(fcs_raw)
gc()

# Save SCE object.
sprintf("Saving SCE object with arcsinh transformed data.")
saveRDS(sce, file.path(results_folder, "sce_arcsinh.rds"))

sprintf("Done")

