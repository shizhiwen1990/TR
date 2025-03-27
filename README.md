data/ – Folder for Independent training set and independent validation set.
DEGs_plot.R – Script for differential expression gene (DEG) visualization.
Deep_learning_for_identifying_TR_CD8+TILs.ipynb – Jupyter Notebook implementing deep learning for identifying tumor-reactive CD8+ TILs.
Deep_learning_for_identifying_TR_CD8+TILs_Revision.ipynb – Jupyter Notebook implementing deep learning for identifying tumor-reactive CD8+ TILs,including QC processes such as confusion matrix and cross validation.
Heatmap_for_cell_type.R – Generates heatmaps to visualize cell type-specific gene expression.
SCTtransform+z_score.R – Performs SCTransform normalization and computes z-scores.
Skyscraper_plots_for_scTCR.R – Creates skyscraper plots for single-cell TCR analysis.
cellchat.R – Runs CellChat for cell-cell communication analysis.
pySCENIC.R – Executes the pySCENIC pipeline for gene regulatory network inference.
pySCENIC_DE_Regulon.R – Identifies differentially expressed regulons using pySCENIC.
scRNA_analysis.R – General scRNA-seq analysis workflow.
scTCR_analysis.R – Pipeline for analyzing single-cell TCR sequencing data.

Prerequisites
R (≥4.0)
Python (if running deep learning models or pySCENIC)
Required R packages: Seurat, CellChat, SCTransform, ggplot2, etc.

License
This project is licensed under the MIT License.
