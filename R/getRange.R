# -------------------------------------------------------------------------------------
# Name             : getRange 
# Description      : Get the next lower and upper limit
# R Version        : 4.1.0 
# -------------------------------------------------------------------------------------
# Author           : Justine B. Dayrit 
# Author Email     : j.bonifacio@irri.org
# Script Version   : 1.0.0
# Modified by      : 
# Date Modified    : 
# Changes made     : 
# -------------------------------------------------------------------------------------
#' @name getRange
#' @aliases getRange
#  @aliases getRange.default
#' @title getRange
#'
#' @description Get the next lower and upper limit
#'
#' @param index reference number
#' @param topToBottom logical, whether plot number will be written from to to bottom
#' @param serpentine logical, whether plot number will be arranged as serpentine order
#' @param numRowPerGrp number of rows per block/replicate
#' @param numColPerGrp number of rows per block/replicate
#' @param numGrpRow number of blocks/replicates per row
#' @param numGrpCol number of blocks/replicates per column
#' @param numRow number of rows per replicate/field
#' @param numCol number of columns per replicate/field
#' @param rowIndexLL previous row lower limit
#' @param rowIndexUL previous row upper limit
#' @param colIndexLL previous column lower limit
#' @param colIndexUL previous column upper limit
#' 
#' @return a set of integers
#' 
# -------------------------------------------------------------------------------------

getRange <- function (index, topToBottom = TRUE, serpentine = FALSE,
                      numRowPerGrp, numColPerGrp, numGrpRow, numGrpCol,numRow, numCol,
                      rowIndexLL, rowIndexUL, colIndexLL, colIndexUL) { #UseMethod("getRange") }

#getRange.default <- function (index, topToBottom = TRUE, serpentine = FALSE,
#                      numRowPerGrp, numColPerGrp, numGrpRow, numGrpCol,numRow, numCol,
#                      rowIndexLL, rowIndexUL, colIndexLL, colIndexUL){
  
  if (topToBottom) {
    
    if(serpentine & numCol != 1){
      
      if (ceiling(index/numGrpRow)%%2 == 0){
        if (index%%numGrpRow != 0) { rowIndex <- index%%numGrpRow 
        } else { rowIndex <- numGrpRow }
        if (rowIndex == 1) { rowIndexUL <- numRow  
        } else { rowIndexUL <- rowIndexUL - numRowPerGrp }
        rowIndexLL <- rowIndexUL - numRowPerGrp + 1
      } else{
        if (index%%numGrpRow != 0) { rowIndex <- index%%numGrpRow 
        } else { rowIndex <- numGrpRow }
        if (rowIndex == 1) { rowIndexLL <- rowIndex  
        } else { rowIndexLL <- rowIndexLL + numRowPerGrp }
        rowIndexUL <- rowIndexLL + numRowPerGrp - 1
      }
      
    } else{
      
      if (index%%numGrpRow != 0) { rowIndex <- index%%numGrpRow 
      } else { rowIndex <- numGrpRow }
      if (rowIndex == 1) { rowIndexLL <- rowIndex  
      } else { rowIndexLL <- rowIndexLL + numRowPerGrp }
      rowIndexUL <- rowIndexLL + numRowPerGrp - 1
      
    }
    
    colIndex <- ceiling(index/numGrpRow)
    colIndexLL <- (colIndex * numColPerGrp) - numColPerGrp + 1
    colIndexUL <- colIndexLL + numColPerGrp - 1
    
  } else {
    
    if(serpentine & numRow != 1){
      
      if (ceiling(index/numGrpCol)%%2 == 0){
        if (index%%numGrpCol != 0) { colIndex <- index%%numGrpCol 
        } else { colIndex <- numGrpCol }
        if (colIndex == 1) { colIndexUL <- numCol } 
        else { colIndexUL <- colIndexUL - numColPerGrp }
        colIndexLL <- colIndexUL - numColPerGrp + 1
      } else{
        if (index%%numGrpCol != 0) { colIndex <- index%%numGrpCol 
        } else { colIndex <- numGrpCol }
        if (colIndex == 1) { colIndexLL <- colIndex  } 
        else { colIndexLL <- colIndexLL + numColPerGrp }
        colIndexUL <- colIndexLL + numColPerGrp - 1
      }
      
    } else{
      
      if (index%%numGrpCol != 0) { colIndex <- index%%numGrpCol 
      } else { colIndex <- numGrpCol }
      if (colIndex == 1) { colIndexLL <- colIndex  } 
      else { colIndexLL <- colIndexLL + numColPerGrp }
      colIndexUL <- colIndexLL + numColPerGrp - 1
      
    }
    
    rowIndex <- ceiling(index/numGrpCol)
    rowIndexLL <- (rowIndex * numRowPerGrp) - numRowPerGrp + 1
    rowIndexUL <- rowIndexLL + numRowPerGrp - 1
    
  }
  
  return(invisible(data.frame(rowIndexLL, rowIndexUL, colIndexLL, colIndexUL)))
  
} # end of getRange function
