
# Vignette for single cell preliminary analysis

## The main source

This Rmarkdown is based on https://satijalab.org/seurat/articles/pbmc3k_tutorial
with some changes in graphics and added exploration of using subset() function
of Seurat to find markers according to subsets of clusters.  This function
may enhance the use of Seurat in other ways too


## What we do

`seurat_jieun_sct.Rmd` is the markdown file with about 25 code chunks.
For practical use, one
can add file saving commands for computed files and figures, and skip chunks with unwanted
figures etc.

The chunks cover all elements necessary in single set analysis, normalization, PCA,
dimension reduction, clustering, marker identification with figues or tables for
every step.

That said, Seurat package has much
more, allowing to adapt to specific needs of a project.

I replaced older normalization and
identification of highly variable features (genes) with a single command that
applies improved transformation published by Satija Lab in 2022.

## Next steps

A script that creates input directory from gene/cell matrix.

Workflow with input and output paths as arguments.

Testing the workflow on an input with size typical for new data sets,
with 15 times more cells and
ca. 2.5 more UMIs per cell on the average, which are realistic parameters
in new data sets.  This will require converting GEO data to a form
required by Seurat and more desktop time...

Improving some figures using ggplot2, including visualization of clusters
for large number of cells (reducing file size etc).

