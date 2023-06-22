# -------------------------------------------------------------------------------------
# Name             : randomizeAlphaLatticeWithRepCheck
# Description      : Generate randomization for Alpha Lattice designs with repeated checks
# R Version        : 4.0.3
# -------------------------------------------------------------------------------------
# Author           : Justine B. Bonifacio
# Author Email     : j.bonifacio@irri.org
# Script Version   : 1.0.0
# Modified by      : 
# Date Modified    : 
# Changes made     : 
# -------------------------------------------------------------------------------------
#' @name randomizeAlphaLatticeWithRepCheck
#' @aliases randomizeAlphaLatticeWithRepCheck
#  @aliases randomizeAlphaLatticeWithRepCheck.default
#' @title Randomization for Alpha Lattice designs with Repeated Checks
#'
#' @description Generate randomization.
#'
#' @param numTrmt number of treatment
#' @param numFieldRow number of field rows
#' @param numFieldCol number of field columns
#' @param numRowPerRep number of row per replicate
#' @param numColPerRep number of column per replicate
#' @param numRowPerBlk number of row per block
#' @param numColPerBlk number of column per block
#' @param trmtRepPerRep number of replicates of treatment group per replicate
#' @param trmtGroup number of treatments per group
#' @param trmtList NULL or vector or character containing the levels of the treatment
#' 
#' @return A list containing dataframe and statistical design array.
#'
#' @examples
#'randomizeAlphaLatticeWithRepCheck <- function(numTrmt = 23, numFieldRow = 20, numFieldCol = 4, 
#'                                              numRowPerRep = 20, numColPerRep = 2, 
#'                                              numRowPerBlk = 10, numColPerBlk = 1,
#'                                              trmtRepPerRep = rep(c(1,3,4), c(16,4,3)),
#'                                              trmtGroup = rep(c(1,2), c(16,7)))
#' 
# -------------------------------------------------------------------------------------

randomizeAlphaLatticeWithRepCheck <- function(numTrmt, numFieldRow, numFieldCol, 
                                              numRowPerRep, numColPerRep, 
                                              numRowPerBlk, numColPerBlk,
                                              trmtRepPerRep, trmtGroup, 
                                              trmtList = NULL) { #UseMethod("randomizeAlphaLatticeWithRepCheck") }
  
#randomizeAlphaLatticeWithRepCheck.default <- function(numTrmt, numFieldRow, numFieldCol,
#                                                      numRowPerRep, numColPerRep,
#                                                      numRowPerBlk, numColPerBlk,
#                                                      trmtRepPerRep, trmtGroup, trmtName = NULL) {
  
  sink(tempfile())
  result <- try(temp <-DiGGer::ibDiGGer(numberOfTreatments = numTrmt,
                                        rowsInDesign = numFieldRow,
                                        columnsInDesign = numFieldCol,
                                        rowsInReplicate = numRowPerRep,
                                        columnsInReplicate = numColPerRep,
                                        rowsInBlock = numRowPerBlk,
                                        columnsInBlock = numColPerBlk,
                                        treatRepPerRep = trmtRepPerRep, 
                                        treatGroup = trmtGroup, 
                                        treatName = trmtList),
                silent = TRUE)
  
  result <- try(temp <-run(temp, continue = TRUE),
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

