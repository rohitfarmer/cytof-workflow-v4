#!/usr/bin/env Rscript --vanilla

# PURPOSE: 1. Calculate fractions of cells per-subject per-cell population.
# 2. Carry out test between male and female samples in COVID day0 set.
# 3. Carry out test between COVID day0 and healthy day0 samples.
# 4. Plot box plots similar to the paper for comparison.

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

# Load external libraries.
cat("Loading external libraries.\n")
suppressMessages(library(tidyverse))

# Assigning results folders (make sure that the script1 is executed before.)
results_folder <- file.path("results", analysis_name)
figures_folder <- file.path("figures", analysis_name)

cell_fractions <- function(dat, cellpop_col){
        print("Calculating per sample per cellpoplation fractions of cells.")
        # Calculate total number of cells per sample.
        df_total_count <- dat %>% group_by(sample_id) %>%
                count(name = "total_count") %>%
                ungroup()

        df_cellpop_count <- dat %>% group_by_at(vars("sample_id", all_of(cellpop_col))) %>%
                count(name = "cellpop_count") %>%
                ungroup()


        # Calculate fraction of resp cells.
        df_fractions <- dplyr::inner_join(df_total_count, df_cellpop_count) %>%
                  dplyr::mutate("fraction" = cellpop_count / total_count)
        return(df_fractions)
}

cellfrac_test_mf <- function(dat, df_cell_fractions, cellpop_col, cond){
# Carry out wilcoxon-test on cell fractions per cell population between male and female.
        print(paste0("Carrying out test between male and female in ", cond))
        dat_t <- dat %>% dplyr::select(sample_id, condition, gender) %>%
                dplyr::filter(condition == cond) %>%
                unique() %>%
                dplyr::inner_join(df_cell_fractions) %>%
                droplevels()

        df_t_out <- tibble()
        for(cp in unique(as.character(dat_t[[cellpop_col]]))){
                x <- dplyr::filter(dat_t, .data[[cellpop_col]] == cp & gender == "F") %>%
                        dplyr::pull(fraction)
                y <- dplyr::filter(dat_t, .data[[cellpop_col]] == cp & gender == "M") %>%
                        dplyr::pull(fraction)
                t_res <- wilcox.test(x, y) %>%
                        broom::tidy()
                df_t_out <- bind_rows(df_t_out, tibble("cellpop" = cp, t_res))
        }
        return(df_t_out)
}

cellfrac_test_cond <- function(dat, df_cell_fractions, cellpop_col, cond1, cond2){
# Carry out wilcoxon-test on cell fractions per cell population between covid and healthy.
        print(paste0("Carrying out test between ", cond1," and ", cond2))
        dat_t <- dat %>% dplyr::select(sample_id, condition, gender) %>%
                dplyr::filter(condition == cond1 | condition == cond2) %>%
                unique() %>%
                dplyr::inner_join(df_cell_fractions) %>%
                droplevels()

        df_t_out <- tibble()
        for(cp in unique(as.character(dat_t[[cellpop_col]]))){
                x <- dplyr::filter(dat_t, .data[[cellpop_col]] == cp & condition == cond1) %>%
                        dplyr::pull(fraction)
                y <- dplyr::filter(dat_t, .data[[cellpop_col]] == cp & condition == cond2) %>%
                        dplyr::pull(fraction)
                t_res <- wilcox.test(x, y) %>%
                        broom::tidy()
                df_t_out <- bind_rows(df_t_out, tibble("cellpop" = cp, t_res))
        }
        return(df_t_out)
}

bp_male_vs_female <- function(dat, cellpop){
        print(paste0("Plotting: ", cellpop))
        samp_gend <- dplyr::filter(dat, condition == "COVID_Day0") %>%
        droplevels() %>%
        dplyr::select(sample_id, gender, condition) %>%
        unique()
        frac <- dplyr::filter(df_cell_fractions, merging1 == cellpop)
        plt_d <- inner_join(samp_gend, frac)
        p <- ggplot(plt_d, aes(x=condition, y=fraction, color = gender)) + 
                  geom_boxplot() +
                    geom_point(position=position_jitterdodge()) +
                  labs(y = paste0(cellpop,"\n(fraction of live cells)"),
                  x = "Condition", color = "Gender") +
                  scale_y_continuous(limits = c(0,0.03))
        ggsave(paste0(cellpop,"_boxplot_mf_covid_day0.png"), plot = p, path = figures_folder, width = 7, height = 5)
}

# ANALYSIS
# Load merging1 exported data.
dat <- readRDS(file.path(results_folder, "sce_post_merging1_export.rds"))

dat <-  dat %>% dplyr::mutate(merging1 = as.character(merging1),
                               merging1 = case_when(merging1 == "Cluster_3" ~ "CD11c+CD16+ DC",
                                                    merging1 == "Cluster_8" ~ "Central memory T helper", 
                                                    merging1 == "Cluster_15" ~ "Mature NK cells",
                                                    merging1 == "Cluster_18" ~ "CD57+PD1+Tc exhausted or senescent",
                                                    TRUE ~ merging1))
#saveRDS(dat, file.path(results_folder, "sce_post_merging1_export_extra_annotations.rds"))
cellpop_col <- "merging1"

# Calculate per-subject per-cluster fractions of cells. 
df_cell_fractions <- cell_fractions(dat, cellpop_col)

# Calcualte wilcoxon-test on fractions of cells between male and female within a condition.
df_w_covidd0 <- cellfrac_test_mf(dat, df_cell_fractions, cellpop_col, "COVID_Day0")
df_w_covidd0_fdr <- df_w_covidd0 %>%
        dplyr::mutate("p.value.fdr" = p.adjust(p.value, method = "fdr"))
#write_tsv(df_w_covidd0_fdr, file.path(results_folder, "male_vs_female_w_covidd0.tsv"))

# Calculate Wilcoxon-test on fractons of cells between samples from covid vs healthy at day 0.
df_w_day0 <- cellfrac_test_cond(dat, df_cell_fractions, cellpop_col, "Healthy_Day0", "COVID_Day0")
df_w_day0_fdr <- df_w_day0 %>%
        dplyr::mutate("p.value.fdr" = p.adjust(p.value, method = "fdr"))
#write_tsv(df_w_day0_fdr, file.path(results_folder, "healthyd0_vs_covidd0_w.tsv"))

# Plot box plot for covid vs healthy.
print("Plotting: covid vs healthy.")
samp_gend <- dplyr::filter(dat, condition == "COVID_Day0" | condition == "Healthy_Day0") %>%
        droplevels() %>%
        dplyr::select(sample_id, gender, condition) %>%
        unique()
frac <- dplyr::filter(df_cell_fractions, merging1 == "pDC")
plt_d <- inner_join(samp_gend, frac)
plt_d$condition <- factor(plt_d$condition, levels = c("COVID_Day0", "Healthy_Day0"))
p <- ggplot(plt_d, aes(x=condition, y=fraction, color = gender)) + 
          geom_boxplot() +
            geom_point(position=position_jitterdodge()) +
          #geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) +
          labs(y = "CD123+ pDCs\n(fraction of live cells)",
          x = "Condition", color = "Gender")
#ggsave("boxplot_covid_vs_healthy_day0.png", plot = p, path = figures_folder, width = 7, height = 5)

# Box plots for male vs female in the 5 cell populations. 
# cellpopulations <- c("CD11c+CD16+ DC", "Central memory T helper", 
#                      "Mature NK cells", "CD57+PD1+Tc exhausted or senescent", "CD11c+ myeloid DC")
 
          # for(cellpop in cellpopulations){
          #         bp_male_vs_female(dat, cellpop)
          # }

# Carefully plot just two populations matching the scale etc.

# CD11c+ DC
samp_gend <- dplyr::filter(dat, condition == "COVID_Day0") %>%
droplevels() %>%
dplyr::select(sample_id, gender, condition) %>%
unique()
frac <- dplyr::filter(df_cell_fractions, merging1 == "CD11c+ myeloid DC")
plt_d <- inner_join(samp_gend, frac)
p <- ggplot(plt_d, aes(x=condition, y=fraction, color = gender)) + 
          geom_boxplot() +
            geom_point(position=position_jitterdodge()) +
          labs(y = paste0("CD11c+ DC","\n(fraction of live cells)"),
          x = "Condition", color = "Gender") +
          scale_y_continuous(limits = c(0,0.03)) +
          theme(text = element_text(size = 20)) 
ggsave(paste0("CD11c+ DC","_boxplot_mf_covid_day0.png"), plot = p, path = figures_folder, width = 7, height = 5)

# CD11c+CD16+ DC
frac <- dplyr::filter(df_cell_fractions, merging1 == "CD11c+CD16+ DC")
plt_d <- inner_join(samp_gend, frac)
p <- ggplot(plt_d, aes(x=condition, y=fraction, color = gender)) + 
          geom_boxplot() +
            geom_point(position=position_jitterdodge()) +
          labs(y = paste0("CD11c+CD16+ DC","\n(fraction of live cells)"),
          x = "Condition", color = "Gender") +
          scale_y_continuous(limits = c(0,0.07)) +
          theme(text = element_text(size = 20)) 
ggsave(paste0("CD11c+CD16+ DC","_boxplot_mf_covid_day0.png"), plot = p, path = figures_folder, width = 7, height = 5)


