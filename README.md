# genomicRegionsFromGeneSymbols

genomicRegionsFromGeneSymbols

This package allows users to retrieve the genomic regions of input gene symbols in the GRCh38/hg38 reference assembly

## Installation via Bioconductor

Once this package will be available on Bioconductor, it will be possibile to install it through the following commands.

Start R (version "4.1") and enter:

```{r, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("genomicRegionsFromGeneSymbols")
```

It will be possible to load the package with the following command:

```{r, eval=FALSE}
library("genomicRegionsFromGeneSymbols")
```


## Usage

The usage of genomicRegionsFromGeneSymbols is very easy. The main function `getGenomicRegionsFromGeneSymbols()` reads a list of gene symbols and some flag parameters and returns the list of genomic regions retrieved from Ensembl.

```{r, eval=TRUE}
neuroblastoma_gene_symbols <- c("AK4", "ALDOC", "EGLN1", "FAM162A", "MTFP1", "PDK1", "PGK1")
neuroblastoma_trackName <- "neuroblastomaPrognosticSignature2020"
gene_symbols_list <- neuroblastoma_gene_symbols
trackName <- neuroblastoma_trackName
verboseFlag <- FALSE
outputFileFormat <- "CSV" # BED or CSV
outputFolder <- "../results/"
saveOutputFileFlag <- FALSE
arrangedAnnotations <- getGenomicRegionsFromGeneSymbols(gene_symbols_list, trackName, saveOutputFileFlag, outputFolder, outputFileFormat, verboseFlag)


```

The function will print all the intermediate messages if `verboseFlag` is set to true, and eventually the `arrangedAnnotations` variable will contain the genomic regions of the query gene symbols

## Contacts

This software was developed by [Davide Chicco](https://www.DavideChicco.it), who can be contacted via email at davidechicco(AT)davidechicco.it