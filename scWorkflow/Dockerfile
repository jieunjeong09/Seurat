FROM rocker/r-ver:4.3.1

# Install R packages
RUN R -e "install.packages(c('Seurat', 'dplyr', 'ggplot2', 'rmarkdown'))"

# Set working directory
WORKDIR /workspace
COPY . /workspace
