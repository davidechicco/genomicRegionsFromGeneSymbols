# test_genomicRegionsFromGeneSymbols <- function() {
test_genomicRegionsFromGeneSymbols <- function() {
    checkEquals(getGenomicRegionsFromGeneSymbols("CDKN3", "test_CDKN3", FALSE, "", "CSV", FALSE), data.frame(GPL96_probeset_ID = "209714_s_at", ensembl_ID = "ENSG00000100526", gene_symbol = "CDKN3", genomic_region = "chr14:54396849-54420218", genomic_region_size_basepairs = 23369, coordinates = "GRCh38/hg38"))
}


