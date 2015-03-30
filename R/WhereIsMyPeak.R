
##' WhereIsMyPeak - Annotate GenomicRanges for gene location
##'
##' @title WhereIsMyPeak
##' @param peaks peaks as GRanges object
##' @param txdb txdb annotation data base
##' @param annotationOrder character vector of prioritising genomic locations for annotation
##' @param typeOrder character vector prioritising gene types if a tied match
##' @param mappingTable mappingTable
##' @param namingTable namingTable
##' @importFrom GenomicFeatures cdsBy
##' @importFrom GenomicFeatures intronsByTranscript
##' @importFrom GenomicFeatures threeUTRsByTranscript
##' @importFrom GenomicFeatures fiveUTRsByTranscript
##' @importFrom BiocParallel bplapply
##' @return list of results
##' @export
##' @author Thomas Schwarzl
WhereIsMyPeak <- function(peaks,
                          annotationOrder = c("CDS", "5' UTR", "3' UTR", "Intron"),
                          typeOrder = c( "non_coding",
                                         "known_ncrna",
                                         "rRNA",
                                         "snoRNA",
                                         "snRNA",
                                         "miRNA", 
                                         "lincRNA", 
                                         "Mt_rRNA", 
                                         "Mt_tRNA",
                                         "sense_intronic",
                                         "transcribed_processed_pseudogene",
                                         "transcribed_unitary_pseudogene",   
                                         "transcribed_unprocessed_pseudogene",
                                         "translated_processed_pseudogene",
                                         "translated_unprocessed_pseudogene",
                                         "unitary_pseudogene",
                                         "unprocessed_pseudogene",
                                         "sense_overlapping",
                                         "3prime_overlapping_ncrna",   
                                         "antisense",
                                         "IG_C_gene",
                                         "IG_D_gene",
                                         "IG_J_gene",
                                         "IG_V_gene",
                                         "IG_C_pseudogene", 
                                         "IG_J_pseudogene",
                                         "IG_V_pseudogene", 
                                         "polymorphic_pseudogene",
                                         "misc_RNA",
                                         "TEC",
                                         "TR_C_gene",
                                         "TR_D_gene",
                                         "TR_J_gene",
                                         "TR_J_pseudogene",
                                         "TR_V_gene",
                                         "TR_V_pseudogene",
                                         "pseudogene",
                                         "processed_pseudogene", 
                                         "processed_transcript",
                                         "LRG_gene", 
                                         "protein_coding"),
                          txdb = NULL,
                          mappingTable = NULL,
                          namingTable = NULL) {
   
   # if txdb is empty initate on through biomart query
   txdb <- .initTxdb(txdb)
   
   # initate mappingTable
   mappingTable <- .initMappingTable(mappingTable)
   
   # initate nameingTable
   namingTable <- .initNamingTable(namingTable)
   
   
   
   #.WhereIsMyPeakEnv(txdb)
   # WhereIsMyPeakEnv <- get("WhereIsMyPeak", envir = .GlobalEnv)
   
   
   # the annotation order must be complet
   #stopifnot(all(annotationOrder %in% c("3'UTR", "5'UTR", "Exon", "Intron")))
   #stopifnot(all(annotationOrder %in% c("protein_coding", "ncRNA", "Bla", "Intron")))
   
   # initiate the empty data structure to store all information
   anno <- list( 
      #types = setNames(rep(NA, 4), c("CDS", "Intron", "5' UTR", "3' UTR")),
      info = list( dataset = "hsapiens_gene_ensembl" ,
                   biomart = "ensembl"),
      final = rep("NA", length(peaks))
   )
   
   # initate the regions which should be annotated
   regions <- .initRegions(txdb)
   
   
   # This is the most time consuming step
   anno[["types"]] <- bplapply(regions, function(x, peaks. = peaks) { 
      annotatedPeaksForRegion(peaks., region = x)
   })
   
   txtoid <- suppressWarnings(select(txdb, columns=c("TXNAME"), keys=keys(txdb), select="all", keytype=c("GENEID")))
   
   rownames(txtoid) <- txtoid[,2]
   
   # Decide which annotation is the best
   decision <- lapply(1:length(peaks), function(i) {
      for(an in annotationOrder) {
         # check if annotation is not empty
         if( !is.null(anno$types[[an]]) && !is.na(anno$types[[an]]) ) {
            # if we have a hit
            if( !is.null(anno$types[[an]][[as.character(i)]]) ) {
               hits <- anno$types[[an]][[as.character(i)]]
               gene.ids <- unique(txtoid[hits,1])
               # when more one gene ids are returned
               pri.id  <- list(id = gene.ids[1], type = mappingTable[gene.ids[1],2] )
               .getPriorityID( gene.ids, mappingTable, typeOrder ) 
               return(list(region = an, gene.id = pri.id$id, type = pri.id$type))
            }
         }
      }
   })
   
   # get the final annotation
   anno[["finalDesc"]] <- t(sapply(decision, function(x) {
      if(is.null(x)) { 
         return(rep(NA, ncol(namingTable)))
      } else { 
         apply(namingTable[x[["gene.id"]],], 2, function(y) paste(y, collapse=",") )
      }
   }))
   
   # get the final geneID
   anno[["finalID"]] <- unlist(lapply(decision, function(x) {
      if(is.null(x)) { 
         return("NA")
      } else { 
         x[["gene.id"]]  
      }
   }))
   
   # get the final geneID
   anno[["finalType"]] <- unlist(lapply(decision, function(x) {
      if(is.null(x)) { 
         return("NA")
      } else { 
         x[["type"]]  
      }
   }))
   
   # get the final region
   anno[["finalRegion"]] <- unlist(lapply(decision, function(x) {
      if(is.null(x)) { 
         return("NA")
      } else { 
         x[["region"]]  
      }
   }))
   
   anno[["finalPlotType"]] <- getTypeForPlot(anno[["finalRegion"]], anno[["finalType"]])
   
   annoFrame <- cbind(as.data.frame(peaks), do.call(cbind, anno[c("finalID", "finalRegion", "finalDesc", "finalType", "finalPlotType")]))
   
   return(annoFrame)
}


.getPriorityID <- function( gene.ids, mappingTable, typeOrder ) {
   if( length(gene.ids) > 1 ) { 
      selection <- mappingTable[gene.ids,] 
      ty <- .getPriorityType(selection[,2], typeOrder)
      return(list(id = selection[selection[,2] == ty,1][1], type = ty))
   } else { 
      return(list(id = gene.ids, type = paste(mappingTable[gene.ids,2],collapse = ",")))
   }
}

.getPriorityType <- function( types, typeOrder ) {
   for( ty in typeOrder ) { 
      if( ty %in% types ) {
         return(ty)
      } 
   } 
}

##' @importFrom GenomicFeatures makeTranscriptDbFromBiomart
.initTxdb <- function(txdb) {
   if(is.null(txdb)) {
      txdb <- makeTranscriptDbFromBiomart( biomart         = "ensembl",
                                           dataset         = "hsapiens_gene_ensembl",
                                           transcript_ids  = NULL,
                                           filters         = "",
                                           id_prefix       = "ensembl_",
                                           host            = "www.biomart.org",
                                           port            = 80,
                                           miRBaseBuild    = NA )
   } else { 
      return(txdb)
   }
} 

.initMappingTable <- function(mappingTable) {
   if(is.null(mappingTable)) {
      ensembl <- useDataset("hsapiens_gene_ensembl", mart = useMart("ensembl"))   
      mappingTable <- getBM(c("ensembl_gene_id", "gene_biotype") , filters="", values=T, uniqueRows = TRUE, mart=ensembl)  
      rownames(mappingTable) <- mappingTable[,1]
   } else { 
      return(mappingTable)
   }
}

.initNamingTable <- function(namingTable) {
   if(is.null(namingTable)) {
      ensembl <- useDataset("hsapiens_gene_ensembl", mart = useMart("ensembl"))   
      namingTable <- getBM(c("ensembl_gene_id", "external_gene_name", "description") , filters="", values=T, uniqueRows = TRUE, mart=ensembl)  
      rownames(namingTable) <- namingTable[,1]
      namingTable <- namingTable[,-1]
   } else {
      return(namingTable)
   }
}

.initRegions <- function(txdb) {
   # Define the regions to map
   # if ( exists("regions", envir = WhereIsMyPeakEnv, inherits=FALSE) ) {
   #    regions <- get("regions", envir = WhereIsMyPeakEnv)
   #} else {
   regions <- list()
   regions[["CDS"]]      <- cdsBy(txdb,  by = "tx", use.names = T)
   regions[["Intron"]]   <- intronsByTranscript( txdb, use.names = T ) 
   regions[["5' UTR"]]   <- fiveUTRsByTranscript( txdb, use.names = T)
   regions[["3' UTR"]]   <- threeUTRsByTranscript(txdb, use.names = T)
   #regions[["Exon"]]     <- exonsBy(txdb, by = "tx", use.names = T)
   #}  
   return(regions)
}

##' @importFrom IRanges elementLengths
##' @importFrom IRanges findOverlaps
##' @importFrom S4Vectors queryHits
##' @importFrom S4Vectors subjectHits
##' @importMethodsFrom BiocGenerics unlist
annotatedPeaksForRegion <- function(peaks, region) {
   region.gr <- unlist(region)
   region.len <- elementLengths(region)
   
   # find all overlaps
   peaks.overlaps <- findOverlaps(peaks, region.gr, ignore.strand = FALSE)
   
   # return false for all if no overlaps are found
   if (length(peaks.overlaps) == 0) {
      return(list(rep(FALSE, length(peaks))))
   }
   
   peaks.overlaps.query.hits <- queryHits(peaks.overlaps)
   peaks.overlaps.subject.hits <- subjectHits(peaks.overlaps)
   
   hits <- lapply(setNames(unique(peaks.overlaps.query.hits), unique(peaks.overlaps.query.hits)), function(i) {
      overlap.ids <- peaks.overlaps.subject.hits[peaks.overlaps.query.hits == i]
      names(region.gr[overlap.ids])
   })
   hits
}


getTypeForPlot <- function(finalType, finalRegion) {
   stopifnot(length(finalType) == length(finalRegion))
   
   return(unlist(lapply(1:length(finalType), function(i) { 
      .r <- finalRegion[[i]]
      .t <- finalType[[i]]
      
      if(.t %in% c("protein_coding") ) {
         return(.r)
      } else {
         return(.t)
      }
   }) ))
}




