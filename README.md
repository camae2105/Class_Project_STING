# STING Activation Signature Project
Scripts and workflows for RNA-seq processing and downstream analysis to build a STING activation signature.

## Overview
This repository contains the supplementary materials for a Digital Biology: Bioinformatics class project. 
It includes the full analysis pipeline used to process RNA-seq datasets and refine a preliminary STING activation signature.

## Repository Structure
- **Scripts/**: All SLURM job scripts and Python scripts used in the pipeline  
- **Metadata_examples/**: Example metadata files for differential expression analysis  
- **Heatmaps/**: Heatmaps generated for the 155-gene and 37-gene signatures  
- **Heatmaps/<dataset>_zscore.tsv**: Z-score matrices used to generate the heatmaps (included within the Heatmaps folder for each dataset)

The folder containing the 155-gene heatmaps also includes the gene list obtained after cluster-based filtering, which was used to overlap the activation, knockout, and rescue experiments to generate the refined 37-gene signature.
