# -------------------------------------------------------------------------------------
# Name             : writePlotNumInfo 
# Description      : Utility functions for layoutIBD
# R Version        : 3.5.1 
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 1.0.0
# Modified by      : 
# Date Modified    : 
# Changes made     : 
# -------------------------------------------------------------------------------------
#' @name writePlotNumInfo
#' @aliases writePlotNumInfo
#  @aliases writePlotNumInfo.default
#' @title writePlotNumInfo
#'
#' @description Utility function.
#'
#' @param info information to be 
#' @param table a matrix
#' @param tableGrp NULL or a matrix
#' @param numGrpRow a numerical value
#' @param numGrpCol a numerical value
#' @param numRowPerGrp a numerical value
#' @param numColPerGrp a numerical value
#' @param incrementWithin a logical value
#' @param increment NULL or a numerical value
#' @param topToBottom logical, whether plot number will be written from to to bottom
#' 
#' @return a matrix
#' 
# -------------------------------------------------------------------------------------

writePlotNumInfo <- function(info, table, tableGrp = NULL,
                             numGrpRow, numGrpCol, numRowPerGrp, numColPerGrp, 
                             incrementWithin = TRUE, increment = NULL,
                             topToBottom = TRUE) { #UseMethod("writePlotNumInfo") }

#writePlotNumInfo.default <- function(info, table, tableGrp = NULL, 
#                                     numGrpRow, numGrpCol, numRowPerGrp, numColPerGrp, 
#                                     incrementWithin = TRUE, increment = NULL,
#                                     topToBottom = TRUE) {
  
  indexCode <- 1
  
  if (topToBottom) {
    for (i in 1:numGrpCol) {
      for (j in 1:numGrpRow) {
        rowIndexLL <- (j * numRowPerGrp) - numRowPerGrp + 1
        rowIndexUL <- rowIndexLL + numRowPerGrp - 1
        colIndexLL <- (i * numColPerGrp) - numColPerGrp + 1
        colIndexUL <- colIndexLL + numColPerGrp - 1
        
        if (incrementWithin) {
          tableGrp[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- indexCode
          table[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- info
          indexCode <- indexCode + 1
          info <- info + increment
        } else {
          table[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- table[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] + info
        }
        
      } ## end stmt -- for (j in 1:numGrpRow)
    } ## end stmt -- for (i in 1:numGrpCol)
  } else {
   for (i in 1:numGrpRow) {
     for (j in 1:numGrpCol) {
       rowIndexLL <- (i * numRowPerGrp) - numRowPerGrp + 1
       rowIndexUL <- rowIndexLL + numRowPerGrp - 1
       colIndexLL <- (j * numColPerGrp) - numColPerGrp + 1
       colIndexUL <- colIndexLL + numColPerGrp - 1

       if (incrementWithin) {
         tableGrp[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- indexCode
         table[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- info
         indexCode <- indexCode + 1
         info <- info + increment
       } else {
         table[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- table[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] + info
       }

     } ## end stmt -- for (z in 1:numBlkRow)
   } ## end stmt -- for (y in 1:numBlkCol)
  }
  
  return(list(table = table, tableGrp = tableGrp))  
  
} ## end of writePlotNumInfo function
