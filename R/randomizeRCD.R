# -------------------------------------------------------------------------------------
# Name             : randomizeRCD
# Description      : Generate randomization for Row-Column Design
# R Version        : 3.6.3 
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 1.0.1
# Modified by      : Justine B. Bonifacio
# Date Modified    : 2021.01.27
# Changes made     : added sink() to hide process output from DiGGer
# -------------------------------------------------------------------------------------
#' @name randomizeRCD
#' @aliases randomizeRCD
#  @aliases randomizeRCD.default
#' @title Randomization for Row Column design
#'
#' @description Generate randomization.
#'
#' @param numTrmt number of treatment
#' @param numRep number of replicates
#' @param numRowPerRep number of row per replicate
#' @param numColPerRep number of column per replicate
#' @param numFieldRow number of field rows
#' @param numFieldCol number of field columns
#' @param trmtList NULL or vector or character containing the levels of the treatment
#' 
#' @return A list containing dataframe and statistical design array.
#'
#' @examples
#' randomizeRCD(numTrmt = 20, numFieldRow = 16, numFieldCol = 5, 
#'              numRep = 4, numRowPerRep = 4, numColPerRep = 5) 
#' 
# -------------------------------------------------------------------------------------

randomizeRCD <- function(numTrmt, numFieldRow, numFieldCol, 
                         numRep, numRowPerRep, numColPerRep,
                         trmtList = NULL) { #UseMethod("randomizeRCD") }

#randomizeRCD.default <- function(numTrmt, numFieldRow, numFieldCol, 
#                                 numRep, numRowPerRep, numColPerRep,
#                                 trmtList = NULL) {
  
  sink(tempfile())
  result <- try(temp <- DiGGer::rcDiGGer(numberOfTreatments = numTrmt, 
                                         rowsInDesign = numFieldRow, 
                                         columnsInDesign = numFieldCol,
                                         rowsInReplicate = numRowPerRep, 
                                         columnsInReplicate = numColPerRep,
                                         treatName = trmtList),
                silent = TRUE)
  sink()
  
  if (all(class(result) == "try-error")) {
    msg <- trimws(strsplit(result, ":")[[1]])
    msg <- trimws(paste(strsplit(msg, "\n")[[length(msg)]], collapse = ""))
    stop(paste("Error in DiGGer:", msg, sep = ""))
  }
  
  setCorrelation(temp, phasenumber = c(1, numRowPerRep, numColPerRep)) ## --- added by AAGulles c/o VIBartolome 2014.05.27
  capture.output(temp)
  
  if (all(class(result) == "try-error")) {
    msg <- trimws(strsplit(result, ":")[[1]])
    msg <- trimws(paste(strsplit(msg, "\n")[[length(msg)]], collapse = ""))
    stop(paste("Error in DiGGer:", msg, sep = ""))
  }
  
  trmtLayout <- getDesign(temp)  
  fieldbook <- temp$dlist
  
  if (!is.null(trmtList)) { trmtLayout <- matrix(fieldbook$ID, nrow(trmtLayout), ncol(trmtLayout)) } 
  fieldbook <- fieldbook[,c("UNIT", "ID", "ENTRY", "ROW", "RANGE", "REP", "TRT")]
  
  return(invisible(list(fieldbook = fieldbook, plan = trmtLayout)))
  
}

