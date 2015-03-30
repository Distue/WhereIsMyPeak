##' @importFrom AnnotationDbi get
.WhereIsMyPeakEnv <- function(txdb) {
   if ( ! exists("WhereIsMyPeakEnv", envir=.GlobalEnv) ) {
      assign("WhereIsMyPeakEnv", new.env(), .GlobalEnv)
   }
   
   WhereIsMyPeakEnv <- get("WhereIsMyPeakEnv", envir=.GlobalEnv)
   
   if (!exists("TXDB", envir = WhereIsMyPeakEnv, inherits = FALSE)) {
      assign("TXDB", txdb, envir = WhereIsMyPeakEnv)
   } else {
      TXDB <- get("TXDB", envir  = WhereIsMyPeakEnv)
      m1 <- unlist(metadata(TXDB))
      m2 <- unlist(metadata(txdb))
      m1 <- m1[!is.na(m1)]
      m2 <- m2[!is.na(m2)]
      
      if ( length(m1) != length(m2) || any(m1 != m2) ) {
         rm(WhereIsMyPeakEnv)
         assign("WhereIsMyPeakEnv", new.env(), .GlobalEnv)
         WhereIsMyPeakEnv <- get("WhereIsMyPeakEnv", envir=.GlobalEnv)
         assign("TXDB", txdb, envir=WhereIsMyPeakEnv)
      }
   }
}
