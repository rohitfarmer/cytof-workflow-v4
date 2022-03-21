#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 4. This script plots
# figures related to clustering for phenotyping data.
# PUBLICATION: https://f1000research.com/articles/6-748

# Read command line arguments passed to the script. 
cmd_args = commandArgs(trailingOnly=TRUE)

# Read yaml file. 
suppressMessages(library(yaml))

# For interactive mode
yaml_file = "pheno_covid_flu_all_gender_100k.yaml"

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
meta_string <- yam$meta_string

# Load external libraries.
cat("Loading external libraries.\n")
suppressMessages(library(cytofWorkflow))
suppressMessages(library(tidyverse))

# Assigning results folders (make sure that the script1 is executed before.)
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

# Load SCE with clustering results saved in script2.
sprintf("Loading SCE with clustering results saved in script2.")
sce <- readRDS(file.path(results_folder, "sce_clust.rds"))

# Extract marker names from daf. There are multiple places from where it can be
# extracted. I am extracting it from colnames of SOM_codes matrix.
marker_names  <- colnames(sce@metadata$SOM_codes)

# Figure 6. Heatmap of the meadian marker intensities across the 20 cell populations.
# obtained with FlowSOM after the metaclustering step with ConsensusClusterPlus.
sprintf("Generating figure6.")
pdf(file =  file.path(figures_folder, "figure6.pdf"), width = 11, height = 8.5)
plotExprHeatmap(sce, features = "type", by = "cluster_id", k = meta_string, 
                        bars = TRUE, perc = TRUE)
dev.off()

# Figure 7. Distributions of marker intensities (arcsinh-transformed) 
# in the 20 cell populations obtained with FlowSOM after the metaclustering step 
# with ConsensusClusterPlus.
sprintf("Generating figure7.")
# pdf(file =  file.path(figures_folder, "figure7.pdf"), width = 11, height = 8.5)
fig7 <- plotClusterExprs(sce, k = meta_string, features = "type")  
ggsave("figure7.pdf", plot = fig7, device = "pdf", path = figures_folder, width = 11, height = 8.5, units = "in")
# dev.off()


# NOTE: Figure 8. Doesn't make sence for phenotyping data since the intention
# is to plot median marker intensities of all the lineage markers and one state
# maker next to them. 
# Figure 8. Heatmap of the median marker intensities of the 10 lineage markers
# and one signaling marker (pS6) across the 20 cell populations obtained with 
# FlowSOM after the metaclustering step with ConsensusClusterPlus (PBMC data).
# sprintf("Generating figure8.")
# for (i in marker_names){
#         logdebug("Generating plot for %s", i)
#         pdf(file = file.path(figures_folder, paste("figure8_", i, ".pdf", sep = "")), width = 20, height = 8.5)
#         fig8 <- plotClusterHeatmap(daf, hm2 = i, k = meta_string, draw_freqs = TRUE)
#         logdebug("Done plotClusterHeatmap for %s", i)
#         dev.off()
#         logdebug("Done dev off for %s", i)
#         rm(fig8)
#         gc()
# }

sprintf("Done")
