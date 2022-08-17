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

arrangedAnnotations <- getGenomicRegionsFromGeneSymbols(gene_symbols_list, trackName,
    saveOutputFileFlag, outputFolder, outputFileFormat, verboseFlag)
```

The function will print all the intermediate messages if `verboseFlag` is set to true, and eventually the `arrangedAnnotations` variable will contain the genomic regions of the query gene symbols. Expected output:

    Query gene symbols' coordinates in the GRCh38/hg38 reference genome:
    gene_symbol           genomic_region      ensembl_ID
    1       AK4   chr1:65147549-65232145   ENSG00000162433  
    2       EGLN1 chr1:231363751-231422287 ENSG00000135766
    3             chr1:231363797-231528603 ENSG00000287856
    4       ALDOC chr17:28573115-28576948  ENSG00000109107
    6       PDK1  chr2:172555373-172608669 ENSG00000152256
    7     FAM162A chr3:122384161-122412334 ENSG00000114023
    5        PGK1 chrX:77910739-78129295  ENSG00000102144


## Contacts

This software was developed by [Davide Chicco](https://www.DavideChicco.it), who can be contacted via email at davidechicco(AT)davidechicco.it
