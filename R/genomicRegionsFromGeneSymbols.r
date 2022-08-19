#' Function that peforms lookup check from jetset
#'
#' @param query query
#' @param result result
#' @return nothing
#' @examples
#' gene_symbols_list1 <- c("AK4", "ALDOC", "EGLN1", "FAM162A", "MTFP1", "PDK1", "PGK1")
#' gene_symbols_list2 <- c("AK4", "ALDOC", "EGLN1", "FAM162A", "MTFP1", "PDK1", "PGK1")
#' jetsetCheckLookup(gene_symbols_list1, gene_symbols_list2)
jetsetCheckLookup <- function(query, result) {
    noResult <- query[which(is.na(result))]
    if(length(noResult) > 0) {
    all.noResult <- paste(noResult, collapse = ', ')
    warning("these items are not recognized: ", all.noResult, call. = FALSE)
    }
    wh.ambig <- which(vapply(result, FUN=length, FUN.VALUE=integer(1)) > 1)
    if(length(wh.ambig) > 0) {
    merged <- vapply(result[wh.ambig], paste, collapse = ',', integer(1))
    all.ambig <- paste(query[wh.ambig], '={', merged, '}', sep = '', collapse = '; ')
    warning("these items map to multiple Entrez IDs: ", all.ambig, call. = FALSE)
    }
}

#' Function that maps gene symbols to probesets of the hgu133a platform from jmap
#'
#' @param symbols_list list of gene symbols
#' @return a list of probesets related to the input gene symbols
#' @import org.Hs.eg.db
#' @import jetset
#' @import AnnotationDbi
#' @examples
#' gene_symbols_list <- c("AK4", "ALDOC", "EGLN1", "FAM162A", "MTFP1", "PDK1", "PGK1")
#' mappings <- (gene_symbols_list)
jmapHgu133aReimplemented <- function(symbols_list) {
  
  nArgs <- sum( !missing(symbols_list))
  if(nArgs == 0) stop("symbols_list must be specified")
  
   nms <- symbols_list
   eg.list <- AnnotationDbi::mget(symbols_list, org.Hs.eg.db::org.Hs.egSYMBOL2EG, ifnotfound = list(NA))
  
   jetsetCheckLookup(nms, eg.list)
   
   eg <- c()
   symb_count <- 1
   for(symb_count in seq_len(length(eg.list))) {
   
        eg[symb_count] <- ifelse((length(eg.list[symb_count]) == 1), eg.list[symb_count], NA)
   }
   
   # eg <- vapply(eg.list, FUN=function(x) if(length(x) == 1) x else NA, FUN.VALUE=character(1))
   # eg_original <- sapply(eg.list, function(x) if(length(x) == 1) x else NA)
   eg <- eg %>% as.matrix() %>% t()
    
    if (all(is.na(eg))) {
    out <- rep(NA, length(eg))
  } else { 
    scores <- jetset::"scores.hgu133a"
    wh <- match(eg, scores$EntrezID)
    out <- rownames(scores)[wh]
    out[is.na(eg)] <- NA
  }
  noProbeset <- nms[is.na(out) & !is.na(eg)]
  if(length(noProbeset) > 0) {
    warning("no probe sets were found for: ", paste(noProbeset, collapse = ", "), call. = FALSE)
  }
  names(out) <- nms
  return(out)
  
  }

#' Function that reads in a list of probesets' ID's for the GPL96 platform and returns the annotations for them, if found
#'
#' @param thisGeneSymbolList list of gene symbols
#' @param thisTrackName name of the track to print in the file
#' @param saveOutputFileFlag says if the output file should be saved or not
#' @param outputFolderName name of the output file directory
#' @param outputFileFormat parameter that specifies the file format of the output
#' @param verboseFlag if TRUE, intermediate messages are going to be printed
#' @export 
#' @import dplyr 
#' @import biomaRt
#' @import jetset
#' @import xml2
#' @import geneExpressionFromGEO
#' @return a dataframe with the annotations containing the chromosome regions
#' @examples
#' gene_symbols_list <- c("AK4", "ALDOC", "EGLN1", "FAM162A", "MTFP1", "PDK1", "PGK1")
#' verboseFlag <- FALSE
#' trackName <- "neuroblastomaPrognosticSignature2020"
#' outputFolder  <- "./test"
#' outputFileFormat <- "BED"
#' outputFileFlag <- FALSE
#' arrangedAnnotations <- getGenomicRegionsFromGeneSymbols(gene_symbols_list, trackName, outputFileFlag, outputFolder, outputFileFormat, verboseFlag)
getGenomicRegionsFromGeneSymbols <- function(thisGeneSymbolList, thisTrackName, saveOutputFileFlag,  outputFolderName, outputFileFormat, verboseFlag) {
  
  if(saveOutputFileFlag == TRUE & outputFileFormat != "CSV" & outputFileFormat != "BED") {
  
    cat("Error: the outputFileFormat should be CSV or BED, while it is ", outputFileFormat, "\n", sep="")
    cat("The program will stop here\n")
    return(NULL)
  }
  
  # For the printed files
  num_to_return <- 1
  exe_num <-  sample(seq_len(as.numeric(10000)), num_to_return)
  
  outputFileName <- ""
  if(saveOutputFileFlag) dir.create(outputFolderName)
  if(verboseFlag & saveOutputFileFlag) cat("created folder ", outputFolderName, "\n", sep="")
  if(saveOutputFileFlag) outputFileName <- paste0(outputFolderName, "/", thisTrackName, "_rand", exe_num,".", tolower(outputFileFormat))
  
  flank_len <- 25
  flank_stream <- "downstream"
  flank_type <- "gene"
  
  # jmap_output <- jetset::jmap("hgu133a", symbol = thisGeneSymbolList))
  jmap_output <- jmapHgu133aReimplemented(thisGeneSymbolList)
  
  probesets_jmap_original  <- jmap_output %>% as.data.frame() %>% as.list()
  probesets_jmap_original <- probesets_jmap_original$.

  genesWithoutProbesets <- thisGeneSymbolList[probesets_jmap_original %>% is.na() %>% which()]
  if(verboseFlag) { 
    if(genesWithoutProbesets %>% length() >= 1) {
    cat("probesets not found for the following genes:\n")
    print(genesWithoutProbesets)
    } else {
    cat("we retrieved a GPL96 probeset for each of the ", thisGeneSymbolList %>% length(), " genes\n", sep="")
    }
  }

  probeset_array  <- probesets_jmap_original[!(probesets_jmap_original %>% is.na())]

  if(verboseFlag) cat("GPL96 probesets (", probeset_array %>% length()," probesets found): \t", sep="")
  if(verboseFlag) print(probeset_array)

  stopifnot(flank_stream %in% c("downstream", "upstream"))
  stopifnot(flank_type %in% c("gene", "coding_gene", "transcript", "coding_transcript"))
  flank_stream <- paste0(flank_stream, "_flank")
  flank_type <- paste0(flank_type, "_flank")
  
  # if(verboseFlag) cat("Attention: the function retrieveAnnotationsGPL96() works only for datasets generated with the GPL96 platform\n")
  # if(verboseFlag) cat("GPL96: [HG-U133A] Affymetrix Human Genome U133A Array ( https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL96 )\n")    
  
  checked_html_text <- "EMPTY_STRING"
  ensembl_URL <- "https://www.ensembl.org/?redirect=no"
  checked_html_text <- xml2::read_html(ensembl_URL)
  checked_url <- "EMPTY_STRING"
  checked_url <- geneExpressionFromGEO::readUrl(ensembl_URL)
  
    if(all(checked_html_text == "EMPTY_STRING"))   {
         
                    cat("The web url ", ensembl_URL," is unavailable right now. Please try again later. The function will stop here\n")
                    return(NULL)
                    
    } else if(all(checked_url == "EMPTY_STRING" | is.null(checked_url[[1]]) )) {
         
                    cat("The web url ", ensembl_URL," is unavailable right now (Error 404 webpage not found). The GEO code might be wrong. The function will stop here\n", sep="")
           return(NULL)        
                    
    } else {
  
        currSpecieMart <- biomaRt::useMart(biomart = "ensembl",  dataset = 'hsapiens_gene_ensembl')
            
        thisAnnotLookup <- biomaRt::getBM(mart = currSpecieMart, 
            attributes = c( "affy_hg_u133a", "ensembl_gene_id", "gene_biotype", "external_gene_name", "chromosome_name", "start_position", "end_position", flank_type), 
            filter = c("affy_hg_u133a", flank_stream), 
            values = list(probeset_array, flank_len), 
            uniqueRows = TRUE, 
            checkFilters = FALSE, bmHeader = TRUE)
            
        if(verboseFlag) cat("~ : ~ : ~  The coordinates here are for the GRCh38/hg38 reference genome ~ : ~ : ~ \n")
        if(verboseFlag) cat("~ : ~ : ~  More information on GRCh38/hg38: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40 ~ : ~ : ~ \n")
            
        theseAnnotations <- thisAnnotLookup

        theseAnnotations <- theseAnnotations %>% dplyr::relocate("Gene end (bp)", .after = "Gene start (bp)")
        theseAnnotations$"Flank (Gene)" <- NULL
        theseAnnotations$"Gene type" <- NULL

        colnames(theseAnnotations)[which(colnames(theseAnnotations)=="Gene name")] <- "gene_symbol"
        colnames(theseAnnotations)[which(colnames(theseAnnotations)=="Gene stable ID")] <- "ensembl_ID"
        colnames(theseAnnotations)[which(colnames(theseAnnotations)=="Chromosome/scaffold name")] <- "chromosome"
        colnames(theseAnnotations)[which(colnames(theseAnnotations)=="Gene start (bp)")] <- "genomic_region_start_basepairs"
        colnames(theseAnnotations)[which(colnames(theseAnnotations)=="Gene end (bp)")] <- "genomic_region_end_basepairs"
        colnames(theseAnnotations)[which(colnames(theseAnnotations)=="AFFY HG U133A probe")] <- "GPL96_probeset_ID"

        theseAnnotations$"chromosome" <- paste0("chr", theseAnnotations$"chromosome")
        
        theseAnnotations$"genomic_region" <- paste0(theseAnnotations$"chromosome", ":",  theseAnnotations$"genomic_region_start_basepairs", "-", theseAnnotations$"genomic_region_end_basepairs")
        
        theseAnnotations <- theseAnnotations[order(theseAnnotations$"chromosome", theseAnnotations$"genomic_region_start_basepairs", theseAnnotations$"genomic_region_end_basepairs"),]
        
        if(outputFileFormat == "CSV") {
        
            theseAnnotations$"genomic_region_size_basepairs" <- theseAnnotations$"genomic_region_end_basepairs" - theseAnnotations$"genomic_region_start_basepairs"
            
            theseAnnotations$"chromosome" <- NULL
            theseAnnotations$"genomic_region_start_basepairs" <- NULL
            theseAnnotations$"genomic_region_end_basepairs" <- NULL
            theseAnnotations$"coordinates" <- "GRCh38/hg38"
        }  
        
        if(verboseFlag) cat("theseAnnotations:\n")
        if(verboseFlag) print(theseAnnotations)
        if(verboseFlag) cat("~ : ~ : ~  The coordinates here are for the GRCh38/hg38 reference genome ~ : ~ : ~ \n")
        if(verboseFlag) cat("~ : ~ : ~  More information on GRCh38/hg38: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40 ~ : ~ : ~ \n")

        # entrezGeneIDlist <- gconvert(theseAnnotations$"Ensembl Gene ID",  organism = "hsapiens", target="ENTREZGENE_ACC")

        number_genes_with_empty_symbol <-  theseAnnotations$"gene_symbol"[theseAnnotations$"gene_symbol" == ""] %>% length()
        number_genes_with_symbol <-  theseAnnotations$"gene_symbol" %>% unique() %>% length() - number_genes_with_empty_symbol
        
            if(number_genes_with_symbol != thisGeneSymbolList %>% unique() %>% length()) {    
                missing_gene_symbols <- setdiff(thisGeneSymbolList, theseAnnotations$"gene_symbol")
                number_missing_gene_symbols <- missing_gene_symbols %>% length()
                
                if(verboseFlag) cat("Attention: only ", number_genes_with_symbol, " gene symbols out of ", thisGeneSymbolList %>% unique() %>% length(), " have annotations here\n", sep="")

                if(number_missing_gene_symbols == 1 && verboseFlag) cat("Attention: ",number_missing_gene_symbols, " gene symbol does not have annotations in GPL96:\n", sep="")
                else if(number_missing_gene_symbols >= 2 && verboseFlag) cat("Attention ",number_missing_gene_symbols, " gene symbols do not have annotations in GPL96:\n", sep="")
                
                if(verboseFlag) print(missing_gene_symbols)
            }
            
            if(theseAnnotations %>% nrow() >=1) {
                cat("\nQuery gene symbols' coordinates in the GRCh38/hg38 reference genome:\n")
                print(theseAnnotations[,c("gene_symbol", "genomic_region", "ensembl_ID")])
                cat("\n")
            } else {
                cat("No genomic regions found for the input gene symbols\n")
            }
            
            SAVE_FILE <- saveOutputFileFlag
            if(SAVE_FILE == TRUE) {
            if(outputFileFormat == "CSV") {
            write.csv(theseAnnotations, file=outputFileName, row.names=FALSE)
            if(verboseFlag) cat("saved annotations in the file ", outputFileName, "\n")
            } else if(outputFileFormat == "BED")  {
            bedHeader <- paste0("track name=\"",thisTrackName,"\" description=\"GRCh38/hg38 coordinates\"")
            write(bedHeader, file=outputFileName, append=TRUE)
                        
            write.table(theseAnnotations[,c("chromosome", "genomic_region_start_basepairs", "genomic_region_end_basepairs")], file=outputFileName, row.names=FALSE, col.names = FALSE, sep="\t", quote = FALSE, append=TRUE)
            
            if(verboseFlag) cat("saved annotations in the file ", outputFileName, "\n")
            }
            }
            return(theseAnnotations)
    }
}


