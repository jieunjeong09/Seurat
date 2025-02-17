
# Vignette for single cell preliminary analysis

## The main source

This Rmarkdown is based on https://satijalab.org/seurat/articles/pbmc3k_tutorial
with some changes in graphics and added exploration of using subset() function
of Seurat to find markers according to subsets of clusters.  This function
may enhance the use of Seurat in other ways too


## The form

scWorkflow -> Rmarkdown -> seurat_jieun.Rmd is the markdown file with about 25 code chunks.
For practical use, one
can add file saving commands for computed files and figures, and skip chunks with unwanted
figures etc.

The chunks cover all elements necessary in single set analysis, normalization, PCA,
dimension reduction, clustering, marker identification with figues or tables for
every step.

That said, Seurat package has much
more, allowing to adapt to specific needs of a project.

## The next step

Although the source vignette was compiled in 2023, the newer data set have deeper
coverage that require a new methodology like one that authored worked out in 2022.
While this vignette is usable, it does not take full advantage of the increased
UNI counts per cell, from ca. 1000 or less to 2500 or more.  So the next step is
to present workflow modified with new methodology, and one improved visualization.
