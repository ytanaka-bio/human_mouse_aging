# human_mouse_aging
Scripts used in Jeong et al. paper

Section 1- Preprocessing and Seurat analysis-
In the preprocessing section, we provided script code to process single-cell RNA-Seq (scRNA-Seq) FASTQ data and perform Seurat preprocessing steps (seurat_pipeline.R). 
Mouse gene symbols were converted to human orthologs using orthologous mapping using Biomart, followed by integration into a unified Seurat object (combine.R).

Section 2- CellChat analysis-
R code used for performing cellchat analysis is provided in
(A) Human folder- for the young human dataset (human_young.R) and the aged human (cell_chat_human_old.R) dataset
(B) Mouse folder- for young mouse (cell_chat_mouse_young.R) and aged mouse (cell_chat_mouse_old.R) dataset
(C) Combined- combined analysis of human young and aged population (combined_analysis.R) and for mouse young and aged population (mouse_cellchat_combined_analysis.R)

Section 3- Subtype analysis-
Further cell subtype analysis was performed on both human and mouse. Scripts for preprocessing and normalization of cell subtypes are described in Subtype.R file. Further analysis was performed on major cell types such as Oligodendrocyte Precursor cells (OPC), Oligodendrocyte cells (ODCs), Astrocyte, and Microglia (subtype_analysis.R). Monocle analysis (monocle.R) was carried out to infer cellular trajectories in ODCs. Additionally, a separate analysis was performed to clarify the relationship between OPCs and ODCs.
