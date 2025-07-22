# app.R

options(shiny.timeout = 60 * 60)
options(shiny.maxRequestSize = 50 * 1024^2)

# Load Libraries
library(shiny)
library(shinydashboard)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(plotly)
library(AnnotationDbi)
library(matrixStats)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(viridis)
library(DT)
library(biomaRt)

#--------------------------------------#
#                 UI                  #
#--------------------------------------#

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

#--------------------------------------#
#               SERVER                #
#--------------------------------------#

server <- function(input, output, session) {
  
  # Setup Ensembl biomaRt for Mus musculus
  get_mart <- reactive({
    useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  })
  
  # Annotate ENSEMBL transcript IDs with gene symbols
  annotate_genes_biomart <- function(ids, mart) {
    ids_clean <- gsub("\\..*", "", ids)
    mapping <- getBM(
      attributes = c("ensembl_transcript_id", "mgi_symbol"),
      filters = "ensembl_transcript_id",
      values = ids_clean,
      mart = mart
    )
    setNames(mapping$mgi_symbol, mapping$ensembl_transcript_id)
  }
  
  # Load input files
  count_data <- reactive({
    req(input$count_file)
    counts <- read.csv(input$count_file$datapath, row.names = 1)
    counts[] <- round(counts)
    counts
  })
  
  metadata <- reactive({
    req(input$metadata_file)
    read.csv(input$metadata_file$datapath)
  })
  
  # Run DESeq2
  run_deseq2 <- reactive({
    dds <- DESeqDataSetFromMatrix(countData = count_data(), colData = metadata(), design = ~ Condition)
    dds <- DESeq(dds)
    res <- results(dds)
    list(dds = dds, res = res)
  })
  
  # Create annotated DESeq2 table
  render_deseq_table <- function(res_df, mart) {
    res_df$Gene <- gsub("\\..*", "", rownames(res_df))
    res_df$GeneSymbol <- annotate_genes_biomart(rownames(res_df), mart)[res_df$Gene]
    res_df
  }
  
  output$differential_expression_table <- renderDT({
    res_df <- as.data.frame(run_deseq2()$res)
    res_df <- render_deseq_table(res_df, get_mart())
    datatable(res_df, options = list(pageLength = 10, scrollX = TRUE), rownames = TRUE)
  })
  
  output$download_table <- downloadHandler(
    filename = function() paste("DESeq2_Results_", Sys.Date(), ".csv", sep = ""),
    content = function(file) {
      res_df <- as.data.frame(run_deseq2()$res)
      res_df <- render_deseq_table(res_df, get_mart())
      write.csv(res_df, file, row.names = FALSE)
    }
  )
  
  # PCA Plot
  output$pca_plot_interactive <- renderPlotly({
    vsd <- vst(run_deseq2()$dds, blind = FALSE)
    pca <- prcomp(t(assay(vsd)))
    pca_data <- as.data.frame(pca$x)
    pca_data$Condition <- colData(run_deseq2()$dds)$Condition
    p <- ggplot(pca_data, aes(PC1, PC2, color = Condition)) +
      geom_point(size = 3) +
      labs(title = "PCA Plot") +
      theme_minimal()
    ggplotly(p)
  })
  
  # Volcano Plot
  output$volcano_plot_interactive <- renderPlotly({
    df <- as.data.frame(run_deseq2()$res)
    df$Gene <- gsub("\\..*", "", rownames(df))
    df$GeneSymbol <- annotate_genes_biomart(rownames(df), get_mart())[df$Gene]
    p <- ggplot(df, aes(log2FoldChange, -log10(pvalue), color = pvalue < 0.05, text = GeneSymbol)) +
      geom_point(size = 1.5) +
      labs(title = "Volcano Plot") +
      theme_minimal()
    ggplotly(p, tooltip = "text")
  })
  
  # Heatmap
  output$heatmap_plot <- renderPlotly({
    vsd <- vst(run_deseq2()$dds, blind = FALSE)
    mat <- assay(vsd)
    top_genes <- head(order(rowVars(mat), decreasing = TRUE), 20)
    mat <- mat[top_genes, ]
    df <- as.data.frame(as.table(mat))
    colnames(df) <- c("Gene", "Sample", "Expression")
    df$Gene <- annotate_genes_biomart(rownames(mat), get_mart())[df$Gene]
    p <- ggplot(df, aes(Sample, Gene, fill = Expression)) +
      geom_tile() +
      scale_fill_viridis_c() +
      theme_minimal() +
      labs(title = "Top Variable Genes") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(p)
  })
  
  # KEGG Enrichment
  run_kegg_analysis <- reactive({
    res <- run_deseq2()$res
    degs <- rownames(res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ])
    if (length(degs) == 0) return(NULL)
    degs <- gsub("\\..*", "", degs)
    mart <- get_mart()
    gene_map <- getBM(c("ensembl_transcript_id", "entrezgene_id"),
                      filters = "ensembl_transcript_id",
                      values = degs, mart = mart)
    entrez_ids <- na.omit(gene_map$entrezgene_id)
    if (length(entrez_ids) == 0) return(NULL)
    enrichKEGG(gene = entrez_ids, organism = "mmu")
  })
  
  output$kegg_plot_interactive <- renderPlotly({
    kegg_df <- as.data.frame(run_kegg_analysis())
    if (nrow(kegg_df) == 0) return(NULL)
    p <- ggplot(kegg_df, aes(reorder(Description, Count), Count, fill = -log10(p.adjust))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_viridis_c() +
      labs(title = "KEGG Pathway Enrichment", x = "Pathway", y = "Gene Count") +
      theme_minimal()
    ggplotly(p)
  })
  
  output$kegg_table <- renderDT({
    kegg_res <- run_kegg_analysis()
    
    req(!is.null(kegg_res))
    
    kegg_df <- as.data.frame(kegg_res)
    if (nrow(kegg_df) == 0) {
      return(datatable(data.frame(Message = "No KEGG enrichment results found.")))
    }
    
    datatable(kegg_df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$download_kegg <- downloadHandler(
    filename = function() paste0("KEGG_Enrichment_", Sys.Date(), ".csv"),
    content = function(file) {
      kegg_res <- run_kegg_analysis()
      if (is.null(kegg_res)) {
        write.csv(data.frame(Message = "No enrichment results."), file, row.names = FALSE)
      } else {
        write.csv(as.data.frame(kegg_res), file, row.names = FALSE)
      }
    }
  )
}

# Run the app
shinyApp(ui, server)

