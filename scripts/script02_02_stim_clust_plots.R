#!/usr/bin/env Rscript --vanilla

# PURPOSE: To run CyTOF workflow (Robinson's) version 4. This script plots
# figures related to clustering for stimulation samples. In figure 8 it plots
# per state marker. Script2_1_clust_plots.R plots figure 8 per type maker.

# PUBLICATION: https://f1000research.com/articles/6-748

# Read command line arguments passed to the script. 
cmd_args = commandArgs(trailingOnly=TRUE)

# Temp. If using RStudio.
# setwd("/Users/farmerr2/locus/sandbox/projects/mini")

# Read yaml file. 
suppressMessages(library(yaml))

# For interactive mode
yaml_file = "phospho_h1n1_30k_20_10k.yaml"

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
suppressMessages(library(logging))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(flowCL))
suppressMessages(library(stringr))

# LOGGING
log_file <- file.path("logs", paste(analysis_name, ".log", sep = ""))
basicConfig(level='FINEST')
addHandler(writeToFile, file=log_file, level='DEBUG')

# Assigning results folders (make sure that the script1 is executed before.)
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

# Load daFrame with clustering results saved in script2.
loginfo("Loading daFrame with clustering results saved in script2.")
daf <- readRDS(file.path(results_folder, "daf_clust.rds"))

# Extract state marker names from daf.
marker_names  <- state_markers(daf)

# Figure 6. Heatmap of the meadian marker intensities across the 20 cell populations.
# obtained with FlowSOM after the metaclustering step with ConsensusClusterPlus.
loginfo("Generating figure6.")
pdf(file =  file.path(figures_folder, "figure6.pdf"), width = 11, height = 8.5)
fig6 <- plotClusterHeatmap(daf, hm2 = NULL, k = meta_string, m = NULL, cluster_anno = TRUE, draw_freqs = TRUE) 
# exp_mat <- fig6[[1]]@ht_list$expression@matrix
# exp_mat_cols <- colnames(exp_mat)
# clust_anno <- as.character()
# for(i in 1:dim(exp_mat)[1]){
#         loginfo("Annotating cluster %s", i)
#         mark_int <- as.character()
#         for(mark in exp_mat_cols){
#                 mark_exp <- exp_mat[i,mark]
#                 if(mark_exp == 0){
#                         mark_int <- append(mark_int, paste(mark, "-", sep =""))
#                 }else if(mark_exp > 0 & mark_exp <= 0.3){
#                         mark_int <- append(mark_int, paste(mark, "lo", sep = ""))
#                 }else if(mark_exp > 0.3 & mark_exp <= 0.8){
#                         mark_int <- append(mark_int, paste(mark, "+", sep = ""))
#                 }else if(mark_exp > 0.8){
#                         mark_int <- append(mark_int, paste(mark, "hi", sep = ""))
#                 }
#         }
#         mark_str <- paste(mark_int, sep = "", collapse = "")
#         try(anno_res <- flowCL(mark_str))
#         cell_pop <- anno_res$Cell_Label[[mark_str]][[1]]
#         if(length(cell_pop) == 0){
#                 clust_anno <- append(clust_anno, paste("Cluster", i, sep ="_"))
#         }else {
#                 clust_anno <- append(clust_anno, cell_pop)
#         }
# }
# merging1 <- data.frame(original_cluster = 1:dim(exp_mat)[1], new_cluster = clust_anno)
# write.table(merging1, file.path("meta", paste(analysis_name, "merging1.txt", sep = "_")), sep = "\t", row.names = FALSE)
dev.off()
rm(fig6)
gc()

# Figure 7. Distributions of marker intensities (arcsinh-transformed) 
# in the 20 cell populations obtained with FlowSOM after the metaclustering step 
# with ConsensusClusterPlus.
fig_7 <- FALSE
if(fig_7){
loginfo("Generating figure7.")
# pdf(file =  file.path(figures_folder, "figure7.pdf"), width = 11, height = 8.5)
fig7 <- plotClusterExprs(daf, k = meta_string, markers = "type")  
ggsave("figure7.pdf", plot = fig7, device = "pdf", path = figures_folder, width = 11, height = 8.5, units = "in")
# ggsave("figure7.png", plot = fig7, device = "png", path = figures_folder, width = 11, height = 8.5, units = "in")
# dev.off()
rm(fig7)
gc()
}

# Figure 8. Heatmap of the median marker intensities of all the type markers
# and one state marker across the meta clusters. 
loginfo("Generating figure8.")
for (i in marker_names){
        logdebug("Generating plot for %s", i)
        pdf(file = file.path(figures_folder, paste("figure8_", i, ".pdf", sep = "")), width = 20, height = 8.5)
        fig8 <- plotClusterHeatmap(daf, hm2 = i, k = meta_string, draw_freqs = TRUE)
        logdebug("Done plotClusterHeatmap for %s", i)
        dev.off()
        logdebug("Done dev off for %s", i)
        rm(fig8)
        gc()
}

loginfo("Done")
