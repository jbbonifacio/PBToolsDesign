# -------------------------------------------------------------------------------------
# Name             : randomizeIBD
# Description      : Generate randomization for incomplete block designs
# R Version        : 3.6.3
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 3.0.1
# Modified by      : Justine B. Bonifacio
# Date Modified    : 2021.01.27
# Changes made     : added sink() to hide process output from DiGGer
# -------------------------------------------------------------------------------------
#' @name randomizeIBD
#' @aliases randomizeIBD
#  @aliases randomizeIBD.default
#' @title Randomization for Incomplete Block designs
#'
#' @description Generate randomization.
#'
#' @param numTrmt number of treatment
#' @param numRep number of replicates
#' @param numBlk number of blocks per replicate
#' @param blksize number of plots per block
#' @param numRowPerBlk number of row per block
#' @param numColPerBlk number of column per block
#' @param numRowPerRep number of row per replicate
#' @param numColPerRep number of column per replicate
#' @param numFieldRow number of field rows
#' @param numFieldCol number of field columns
#' @param trmtList NULL or vector or character containing the levels of the treatment
#' 
#' @return A list containing dataframe and statistical design array.
#'
#' @examples
#' randomizeIBD(numTrmt = 16, numRep = 3, numRowPerRep = 16, numColPerRep = 1, 
#'              numBlk = 4, blksize = 4, numRowPerBlk = 4, numColPerBlk = 1, 
#'              numFieldRow = 16, numFieldCol = 3) 
#' 
# -------------------------------------------------------------------------------------

randomizeIBD <- function(numTrmt, numFieldRow, numFieldCol, 
                         numRep, numRowPerRep, numColPerRep,
                         numBlk = NULL, blksize = NULL, numRowPerBlk = NULL, numColPerBlk = NULL,
                         trmtList = NULL) { #UseMethod("randomizeIBD") }

#randomizeIBD.default <- function(numTrmt, numFieldRow, numFieldCol, 
#                                 numRep, numRowPerRep, numColPerRep,
#                                 numBlk = NULL, blksize = NULL, numRowPerBlk = NULL, numColPerBlk = NULL,
#                                 trmtList = NULL) {
  
#  if (is.null(numRowPerBlk) & is.null(numColPerBlk)) { block2DInput <- TRUE  
#  } else { block2DInput <- c(numRowPerBlk, numColPerBlk) }
  
  sink(tempfile())
  result <- try(temp <-ibDiGGer(numberOfTreatments = numTrmt,
                                rowsInDesign = numFieldRow,
                                columnsInDesign = numFieldCol,
                                rowsInReplicate = numRowPerRep,
                                columnsInReplicate = numColPerRep,
                                rowsInBlock = numRowPerBlk,
                                columnsInBlock = numColPerBlk,
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
  # if (!is.na(match("B121", names(fieldbook)))) { fieldbook[,match("B121", names(fieldbook))] <- NULL }
  # if (!is.na(match("B111", names(fieldbook)))) {
  #   # --- Design == AlphaLattice
  #   if (!is.null(numRowPerBlk) & !is.null(numColPerBlk)) { names(fieldbook)[match("B111", names(fieldbook))] <- "BLOCK" }
  #   # --- Design == RowColumn
  #   if (is.null(numRowPerBlk) & is.null(numColPerBlk)) { 
  #     fieldbook[,match("B111", names(fieldbook))] <- NULL 
  #   }
  # } ## end stmt -- if (!is.na(match("B111", names(fieldbook))))
  
  return(invisible(list(fieldbook = fieldbook, plan = trmtLayout)))
  
}

