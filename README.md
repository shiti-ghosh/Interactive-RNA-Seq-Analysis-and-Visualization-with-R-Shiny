# RNAseq-Shiny-Dashboard
This repository contains an interactive Shiny dashboard designed for exploratory analysis of RNA-seq data processed through pseudoalignment-based quantification tools such as Kallisto. The platform is optimized for Mus musculus (mouse) datasets and is capable of performing differential gene expression analysis using DESeq2, visualization like PCA Plot, Volcano Plot, Heatmap, and pathway enrichment, from a count matrix, without requiring full genome alignment or raw BAM/SAM files.

The dashboard enables wet-lab researchers and computational biologists to rapidly gain insights into their transcriptomic datasets. With a clean UI and modular features, this app is particularly useful for academic research, pre-publication analysis, or collaborative biology projects involving multiple conditions or time points.

Key Features:

The application is powered by DESeq2 (R) for differential expression analysis, with visualizations built using ggplot2, plotly, and clusterProfiler. Once the user uploads a transcript- or gene-level count matrix and a corresponding metadata table (CSV format), the app enables them to:

1. View differentially expressed genes in an interactive table
2. Explore PCA-based clustering and sample variability
3. Examine Volcano plots of fold changes vs significance
4. Visualize top variable genes in an interactive heatmap
5. Perform KEGG pathway enrichment analysis based on DE genes
6. These results can be downloaded directly from the interface for further interpretation or reporting.


Using the Dashboard Application:

Input Requirements:

The app expects two input files in CSV format:

1. Count Matrix
   Rows: Transcript or Gene IDs (Ensembl IDs recommended)
   Columns: Sample names
   Values: Estimated counts (from Kallisto or Salmon outputs)
   
3. Metadata Table
   Must contain a column named Condition corresponding to experimental groups (e.g., treated vs control)
   Rows should match the column names of the count matrix
   
   Metadata Format:

   <img width="103" height="200" alt="image" src="https://github.com/user-attachments/assets/817ca127-d585-4724-88c5-06107b58abcb" />


The app calculates the differential expression values for the genes and helps in visualizing them using:

1. PCA Plot
2. Volcano Plot
3. Enrichment Pathways

Note: the app may take a few minutes to run in the background and generate results from DESeq2(R). Please be patient.
