
# Vignette for single cell preliminary analysis

## The main source

This Rmarkdown is based on https://satijalab.org/seurat/articles/pbmc3k_tutorial
with some changes in graphics and added exploration of using subset() function
of Seurat to find markers according to subsets of clusters.  This function
may enhance the use of Seurat in other ways too


## The form

`scWorkflow -> Rmarkdown -> seurat_jieun.Rmd` is the markdown file with about 25 code chunks.
For practical use, one
can add file saving commands for computed files and figures, and skip chunks with unwanted
figures etc.

The chunks cover all elements necessary in single set analysis, normalization, PCA,
dimension reduction, clustering, marker identification with figues or tables for
every step.

That said, Seurat package has much
more, allowing to adapt to specific needs of a project.

## scTransform()

`scWorkflow -> Rmarkdown -> seurat_jieun_sct.Rmd` replaces older normalization and
identification of highly variable features (genes) with a single command that
in turn, applies improved transformation published by Satija Lab in 2022.

## The next step

I will the improved workflows on a data file with 15 times more cells and
ca. 2.5 more UMIs per cell on the average, which are realistic parameters
in new data sets.  This will require converting GEO data to a form
required by Seurat and more desktop time...


