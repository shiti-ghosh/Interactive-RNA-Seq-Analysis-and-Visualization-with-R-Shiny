
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
options(shiny.timeout = 60 * 60)
library(shiny)
library(shinydashboard)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(plotly)
library(AnnotationDbi)
library(matrixStats)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(viridis)
library(DT)
library(biomaRt)

ui <- dashboardPage(
  dashboardHeader(
    title = span("RNA-Seq Insights Dashboard", style = "font-size: 18px; font-weight: bold;"),
    titleWidth = 300
  ),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
      menuItem("Differential Expression Tables", tabName = "de_tables", icon = icon("table")),
      menuItem("PCA Plot", tabName = "pca", icon = icon("chart-bar")),
      menuItem("Volcano Plot", tabName = "volcano", icon = icon("chart-line")),
      menuItem("Heatmap", tabName = "heatmap", icon = icon("fire")),
      menuItem("KEGG Enrichment", tabName = "kegg", icon = icon("dna"))
    )
  ),
  dashboardBody(
    fluidRow(
      box(
        width = 12, status = "info", solidHeader = FALSE, collapsible = FALSE,
        p(
          "This dashboard enables rapid, interactive analysis of RNA-seq data quantified using pseudoalignment tools like ",
          strong("Kallisto"), " or ", strong("Salmon"), ". It is designed specifically for ",
          em("Mus musculus"), " (mouse) datasets and accepts transcript- or gene-level count matrices as input. ",
          "Built on DESeq2, the app performs differential expression analysis and generates visualizations such as PCA plots, volcano plots, heatmaps, and KEGG pathway enrichment — all without requiring raw read alignment or command-line tools. ",
          "By simplifying transcriptomic exploration, the platform accelerates biological insight and facilitates collaboration between computational and experimental teams."
        )
      )
    ),
    tabItems(
      tabItem(tabName = "upload",
              p("Upload a count matrix and sample metadata. These files should be in CSV format and match by sample IDs."),
              fileInput("count_file", "Upload Count Matrix (.csv)", accept = ".csv"),
              fileInput("metadata_file", "Upload Metadata (.csv)", accept = ".csv")
      ),
      tabItem(tabName = "de_tables",
              fluidRow(
                box(title = "Differential Expression Table", status = "primary", solidHeader = TRUE, width = 12,
                    p("This table displays genes with significant differential expression based on your design matrix. Sort, filter, or download the results below."),
                    DTOutput("differential_expression_table"),
                    br(),
                    downloadButton("download_table", "Download DE Table")
                )
              )
      ),
      tabItem(tabName = "pca",
              fluidRow(
                box(title = "PCA Plot", status = "primary", solidHeader = TRUE, width = 12,
                    p("This plot visualizes global sample variation using Principal Component Analysis (PCA). It helps identify clustering or separation between conditions or batches."),
                    plotlyOutput("pca_plot_interactive", height = "600px")
                )
              )
      ),
      tabItem(tabName = "volcano",
              fluidRow(
                box(title = "Volcano Plot", status = "primary", solidHeader = TRUE, width = 12,
                    p("This plot shows log2 fold-change vs. significance. Significant genes are easily spotted by their position in the upper corners."),
                    plotlyOutput("volcano_plot_interactive", height = "600px")
                )
              )
      ),
      tabItem(tabName = "heatmap",
              fluidRow(
                box(title = "Heatmap", status = "primary", solidHeader = TRUE, width = 12,
                    p("A clustered heatmap of the top differentially expressed genes. Rows = genes, Columns = samples."),
                    plotlyOutput("heatmap_plot", height = "700px")
                )
              )
      ),
      tabItem(tabName = "kegg",
              
              
              fluidRow(
                box(title = "KEGG Enrichment Plot", status = "primary", solidHeader = TRUE, width = 12,
                    p("This tab shows pathway-level enrichment using KEGG. The bar plot visualizes enriched pathways, and a table provides full results."),
                    plotlyOutput("kegg_plot_interactive", height = "800px")
                )
              ),
              fluidRow(
                box(title = "KEGG Enrichment Table", status = "primary", solidHeader = TRUE, width = 12,
                    DTOutput("kegg_table"),
                    br(),
                    downloadButton("download_kegg", "Download KEGG Results")
                )
              )
      )
    )
  )
)
