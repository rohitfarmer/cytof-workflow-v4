# CyTOF Workflow Version 4
This repository contains ready to use code for implementing CyTOF workflow V4 from Mark D. Robinson's group. The code is mentioned in their paper [ref](https://f1000research.com/articles/6-748/v4). I have tried to break the pipeline down into logical chunks to make it easier to run it analysis wise. It also makes it easier to troubleshoot errors and resume analysis from a crash point. The majority of the figures are ggplot objects; however, there few figures that are list type. I have included the code to save the figures accordingly. All the figures are saved through a variable; therefore, no figure panel will popup in Rstudio or an Xserver if running remotely. Figure names are according to the figure numbers in the paper. I would recommend running the figure generation scripts interactively. You may have a different number of samples; therefore, you may want to change the figure dimensions accordingly. Some of the figures, especially the ones that generate marker intensity distributions, require a lot of memory. I run them on a cluster with memory ranging from 64Gb to 500Gb. 

In addition to the figures, the code also saves SCE object (that carries all the input data and results) in an R data structure file after each major calculation step. To ensure if the calculation breaks, then you don't have to run the entire pipeline from the beginning. You can resume from where the last SCE is saved. I have also implemented code to extract information from SCE object and store them in TSV files.

For feedback and collaborations please write to me at [rohit.farmer@gmail.com](mailto:rohit.farmer@gmail.com).

# Requirements for the CyTOF Workflow V4
This implementation of the pipeline requires two meta files (args and panel) describing the experimental design; a YAML file that defines the run specific arguments to the scripts and a manually curated cluster merging file.

1. **A panel file.** It's a tab delimited file with a column for all the metals used, associated antigens per metal, and to define marke class; type for lineage and state for functional.
  1. An example phenotyping panel is available at meta/panel/phenotyping_panel_v3.txt
  2. An example stimulation panel availabe at meta/panel/phospho_panel_v3.txt
2. **An args file.** It's a tab delimited file with a column for FCS filenames to be used in a specific analysis, a number given to each sample as an ID, conditions to study, and a patient id. Sample ID is unique to each FCS file   whereas a patient id can be redundant if more than two samples belong to the same subject/patient. 
   * A separate args file need to be created for every distinct analysis.
3. In addition to the panel and args file CHI's implementaion of the workflow requires a **YAML file** that specifies the location of data, panel and args file, and parameters to be passed to the workflow scripts.
   * A separate YAML file need to be created for every distinct analysis.
4. **Cluster Merging 1 file.** After the meta clustering step a manually curated merging file is required that gives a cell population name to each of the meta clusters.

# Setup an Example Project Environment 
## Clone this Repository
`git clone https://github.com/rohitfarmer/cytof-workflow-v4.git <project_name>`

## Setup Working Directory
After cloning the repository:
```
cd <project_name>
sh setup.sh
```
## Prepare YAML file
Prepare a YAML file per analysis in the meta folder.
```
analysis_name: pheno_covid_flu_all_gender_100k
data_type: phenotyping # or stimulation
data_location: covid_flu_100k # Within data folder.
args_file: args_covid_flu_pheno_all_gender.txt # Within meta folder.
panel_file: panel/cytex_phenotyping_panel_36c.txt # Within meta folder.
condition_levels: 
  - Healthy_Day0
  - Healthy_Day7
  - COVID_Day0
  - COVID_Day7
FACS: TRUE # TRUE for Flow data.
cofactor: 150 # 5 for CyTOF, and 150 for Flow data.
no_of_clusters: 20
meta_string: meta20
tsne_no_cells: 1000
umap_no_cells: 1000
merging1_file: pheno_covid_flu_all_gender_100k_merging1.txt
```

# Scripts to Run in Order
All the scripts can be run interactively by specifying the path to the YAML file at the very beginning of the script. Scripts can also be run with Rscript and passing the YAML file as a command line argument. In the workflow below the script are shown to run through Rscript.  

## Data Import, Transformation and Diagnostic Plots (Figures 1 to 5)
```
Rscript --vanilla scripts/scr01_00_arcsinh_transform.R meta/file.yaml
Rscript --vanilla scripts/scr01_01_diagnostic_plots.R meta/file.yaml
```
**Output Data**
* R data structure with arcsinh transformed data: `results/analysis_name/sce_arcsinh.rds`

## Clustering and Meta Clustering (Figures 6 to 8)
```
Rscript --vanilla scripts/scr02_00_clustering.R meta/file.yaml
Rscript --vanilla scripts/scr02_01_clustering_plots.R meta/file.yaml 
```
**Output Data**
* R data structure with added clustering data: `results/analysis_name/sce_clust.rds`

## Dimensionality Reduction (Figures 9 to 14)
```
Rscript --vanilla scripts/scr03_00_dimensionality_reduction.R meta/file.yaml
Rscript --vanilla scripts/scr03_01_dr_plots.R meta/file.yaml
```
**Output Data**
* R data structure with added dimensionality reduction data: `results/analysis_name/sce_dr.rds`

## Cluster Merging 1 (Figures 15 to 18)
```
Rscript --vanilla scripts/scr04_00_merging1.R meta/file.yaml
```
**Output Data**
* R data structure with added merging 1 data: `results/analysis_name/daf_merging1.rds`

## Differential Abundance (DA) Analysis using GLMM and edgeR (Figures 20 to 22)
```
Rscript --vanilla scripts/scr05_00_diff_abundance.R meta/file.yaml
```
**Output Data**
* Figure 20 data: `results/analysis_name/fig20_data.tsv`
* TSV files with p and adjusted p-values from GLMM (Figure 22): `results/analysis_name/p_glmm0.tsv, p_glmm1.tsv, p_glmm2.tsv`
* TSV file with p and adjusted p-values from edger (Figure 22): `results/analysis_name/p_edger.tsv`              
* Summary statistics from GLMM and edge runs (Figure 22): `results/analysis_name/summary_glmm0.tsv, summary_glmm1.tsv, summary_glmm2.tsv, summary_edger.tsv`

## Differential State (DS) Analysis using LMM (Figures 23 to 26)
```
Rscript --vanilla scripts/scr06_00_diff_state.R meta/file.yaml
```
**Output Data**
* TSV files with p and adjusted p-values from LMM (Figure 24 and 26): `results/analysis_name/p_ds_lmm1.tsv, p_ds_lmm2.tsv, p_ds_lmm3.tsv, p_ds_lmm4.tsv`
* TSV files with summary statistics from LMM (Figure 24 and 26): `results/analysis_name/summary_ds_lmm1.tsv, summary_ds_lmm2.tsv, summary_ds_lmm3.tsv, summary_ds_lmm4.tsv`
