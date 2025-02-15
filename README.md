
# Package for single cell preliminary analysis

## Input and output

Input is a directory in 10X Genomics format, barcodes.tsv for column names = cell identifiers, features.tsv.gz for rows, genes with pairs of symbols, ENSEMBL (unique) and HUGO (sometimes missing or not unique).

Output directory has two visualization files, and a table with a row for every cell giving coordinates in PCA vi
sualization

## Usage: 

from command line, it will be `Rscript seurat.R path_to_input_directory`

### Preprocessing package with interactions

It is important  to address the issues of MT genes and low quality cells, removing non-informative genes, and pro
vide some pre-analysis statistics to select parameters for this preprocessing.  E.g. the "normal" fraction of MT
 reads/UMIs depends on assay and biological source of cells.  For now it is done with default values.

### Options for other clustering and visualization methods

PCA is not always the best so options for other methods provided by Seurat will be added.
