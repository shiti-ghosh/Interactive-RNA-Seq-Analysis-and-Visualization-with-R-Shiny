FROM rocker/shiny:4.3.1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libglpk-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libboost-all-dev \
    && apt-get clean

# Install Bioconductor and required R packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
    && R -e "BiocManager::install(c('shiny', 'shinydashboard', 'DESeq2', 'ggplot2', 'pheatmap', 'plotly', 'AnnotationDbi', 'matrixStats', 'clusterProfiler', 'org.Mm.eg.db', 'enrichplot', 'viridis', 'DT', 'biomaRt'), ask=FALSE, update=TRUE)"

# Copy your Shiny app
COPY . /srv/shiny-server/

# Set permissions
RUN chown -R shiny:shiny /srv/shiny-server

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
