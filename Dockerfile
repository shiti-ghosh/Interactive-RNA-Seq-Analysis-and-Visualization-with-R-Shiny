# Base image with R and Shiny Server pre-installed
FROM rocker/shiny:4.3.1

# Install system dependencies required for your R packages
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libglpk-dev \
    libxt-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libjpeg-dev \
    libpng-dev \
    libgit2-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor and CRAN packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c( \
        'shiny', 'shinydashboard', 'ggplot2', 'pheatmap', 'plotly', \
        'AnnotationDbi', 'matrixStats', 'clusterProfiler', 'org.Mm.eg.db', \
        'enrichplot', 'viridis', 'DT', 'biomaRt', 'DESeq2' \
    ), ask=FALSE, update=TRUE)"

# Copy your app.R to the Shiny server's default app directory
COPY app_scripts/app.R /srv/shiny-server/app.R

# Expose port (Shiny runs on 3838 by default)
EXPOSE 3838

# Run the Shiny server
CMD ["/usr/bin/shiny-server"]

