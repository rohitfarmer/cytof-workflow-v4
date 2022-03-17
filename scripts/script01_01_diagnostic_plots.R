#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 4. This script plots
# diagnostic plots.
# PUBLICATION: https://f1000research.com/articles/6-748

# Read command line arguments passed to the script. 
cmd_args = commandArgs(trailingOnly=TRUE)

# For interactive mode
yaml_file = "pheno_covid_flu_all_gender.yaml"

# Read yaml file. 
suppressMessages(library(yaml))
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
data_type  <- yam$data_type

# Load external libraries.
cat("Loading external libraries.\n")
suppressMessages(library(cytofWorkflow))
suppressMessages(library(tidyverse))
suppressMessages(library(Cairo))
suppressMessages(library(ComplexHeatmap))

# Assigning results folders (make sure that the script1 is executed before.)
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

# Load daFrame with arcsinh transformed data.
sprintf("Loading daFrame with arcsinh transformed data.")
sce <- readRDS(file.path(results_folder, "sce_arcsinh.rds"))

# DIAGNOSTIC PLOTS
sprintf("DIAGNOSTIC PLOTS")

# Figure 1. Per-sample smoothed densities of marker expression (arcsinh-transformed).
sprintf("Generating figure1.")
fig1  <- plotExprs(sce, color_by = "condition")
fig1$facet$params$ncol <- 5
ggsave("figure1.png", plot = fig1, device = "png", path = figures_folder, 
       width = 12, height = 10, units = "in", dpi = 300)
rm(fig1)
gc()

# Figure 2. Barplot showing the number of cells measured for each sample in the PBMC dataset.
sprintf("Generating figure2.")
sprintf("Number of cells per sample.")
n_cells(sce)
fig2 <- plotCounts(sce, color_by = "condition")
ggsave("figure2.png", plot = fig2, device = "png", path = figures_folder, 
       width = 12, height = 10, units = "in", dpi = 300)
rm(fig2)
gc()

# Figure 3. MDS plot.
sprintf("Generating figure3.")
fig3 <- CATALYST::plotMDS(sce, color_by = "condition")
ggsave("figure3.png", plot = fig3, device = "png", path = figures_folder, 
       width = 12, height = 12, units = "in", dpi = 300)
rm(fig3)
gc()

# Figure 4. Heatmap of the median (archsinh-transformed) marker expression 
# of lineage markers and functional markers across all cells measured for 
# each sample in the PBMC dataset.
# For phenotyping data it will not produce any figure because there are no
# functional markers. 
if(data_type == "stimulation"){
        sprintf("Stimulation data. Generating figure4.")
        pdf(file =  file.path(figures_folder, "figure4.pdf"), width = 16, height = 8.5)
        fig4 <- plotExprHeatmap(sce, bin_anno = TRUE, row_anno = TRUE)
        draw(fig4)
        dev.off()
        rm(fig4)
        gc()
}else if(data_type == "phenotyping"){
        sprintf("Phenotyping data skipping figure 4 generation.")
}

# Figure 5. Non-redundancy scores for each of the markers and all samples in the PBMC dataset.
sprintf("Generating figure5.")
fig5 <- plotNRS(sce, features = type_markers(sce), color_by = "condition")
ggsave("figure5.pdf", plot = fig5, device = "pdf", path = figures_folder, width = 12, height = 10, units = "in")
rm(fig5)
gc()

loginfo("Done")
