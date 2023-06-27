# -------------------------------------------------------------------------------------
# Name             : randomizeAlphaLatticeEBS
# Description      : Generate randomization for Alpha Lattice design (from ebsRtools)
# R Version        : 4.0.1
# -------------------------------------------------------------------------------------
# Author           : Justine B. Bonifacio 
# Author Email     : j.bonifacio@irri.org
# Script Version   : 1.1.0
# Modified by      : Justine B. Bonifacio
# Date Modified    : 
# Changes made     : 
# -------------------------------------------------------------------------------------
#' @name randomizeAlphaLatticeEBS
#' @aliases randomizeAlphaLatticeEBS
#  @aliases randomizeAlphaLatticeEBS.default
#' @title Randomization for Alpha Lattice design using ebsRTools package
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
#' @param numBlkRow number of blocks per row
#' @param numBlkCol number of blocks per column 
#' @param numRepRow number of replicates per row
#' @param numRepCol number of replicates per column
#' @param numFieldRow number of field rows
#' @param numFieldCol number of field columns
#' @param topToBottom a logical variable indicating whether plot number will be written from to to bottom
#' @param serpentine a logical variable indicating whether plot number will be arranged as serpentine order
#' @param trmtList NULL or vector or character containing the levels of the treatment
#' @param genLayout logical, whether a layout of the design will be generated
#' 
#' @return A list containing dataframe and statistical design array.
#'
#' @examples
#' randomizeAlphaLatticeEBS (numTrmt = 24, numFieldRow = 8, numFieldCol = 6, 
#'                           numRep = 2, numRowPerRep = 4, numColPerRep = 6, 
#'                           numBlk = 6, blksize = 4, numRowPerBlk = 2, numColPerBlk = 2,
#'                           numBlkRow = 2, numBlkCol = 3, numRepRow = 2, numRepCol = 1,
#'                           topToBottom = TRUE, serpentine = TRUE, trmtList = paste("Entry",1:24,sep=""),
#'                           genLayout = FALSE)
#' 
# -------------------------------------------------------------------------------------

randomizeAlphaLatticeEBS <- function(numTrmt, numFieldRow, numFieldCol,
                                     numRep, numRowPerRep, numColPerRep,
                                     numBlk, blksize, numRowPerBlk, numColPerBlk,
                                     numBlkRow, numBlkCol, numRepRow, numRepCol,
                                     topToBottom = NULL, serpentine = NULL, trmtList = NULL,
                                     genLayout = FALSE) {
  
  sink(tempfile())
  result <- try(temp <- suppressWarnings(randALPHA(trt = 1:numTrmt,
                                                   k = blksize,
                                                   r = numRep+1,
                                                   tag = nchar(numTrmt))),
                silent = TRUE)
  sink()
  
  if (all(class(result) == "try-error")) {
    stop(paste("Error found in using randALPHA from ebsRTools package"))
  }
  
  designparam <- result[[1]]$parameters$design
  
  fieldbookLayout <- result[[1]]$book[1:(numTrmt*numRep),1:3]
  fieldbookEntry <- result[[1]]$book[(1+numTrmt):nrow(result[[1]]$book),4]
  fieldbook <- cbind(fieldbookLayout,entry_id = fieldbookEntry)
  
  #reassign block number from 1 to number of blocks per rep
  fieldbook$block <- rep(rep(1:numBlk,each = blksize),numRep)
  
  #rearrange fieldbook column order
  fieldbook <- fieldbook[,c("plot_number", "block", "entry_id", "replicate")]
  
  #rename fieldbook columns
  names(fieldbook) <- c("PlotNumber","Block","ENTRY","REP")
  
  #merge treatment list with fielbook by entry number
  trmtList <- data.frame("ENTRY" = 1:numTrmt, "ID" = trmtList)
  fieldbook <- merge(fieldbook, trmtList, all.x = TRUE)
  fieldbook <- fieldbook[order(as.numeric(fieldbook$PlotNumber)),]
  rownames(fieldbook) <- 1:nrow(fieldbook)
  
  if(!genLayout){
    tempLayout <- data.frame(index = rep(1:numTrmt,numRep),fieldbook)
    tempLayoutTrmt <- reshape(tempLayout[,c(1,5:6)], v.names = "ID", idvar = "index",
                              timevar = "REP", direction = "wide")
    tempLayoutTrmt <- tempLayoutTrmt[2:ncol(tempLayoutTrmt)]
    names(tempLayoutTrmt) <- paste("Rep",1:numRep, sep = "")
    
    tempLayoutBlock <- reshape(tempLayout[,c(1,4:5)], v.names = "Block", idvar = "index",
                               timevar = "REP", direction = "wide")
    tempLayoutBlock <- tempLayoutBlock[2:ncol(tempLayoutBlock)]
    names(tempLayoutBlock) <- paste("Rep",1:numRep, sep = "")
    
    trmtLayout <- list(tempLayoutTrmt, tempLayoutBlock)
    
    fieldbook <- fieldbook[,c("ID", "ENTRY", "REP", "PlotNumber", "Block")]
    
  } else {
    fieldbook$REP <- as.numeric(fieldbook$REP)
    
    tempLayoutTrmt <- matrix(0, nrow = numFieldRow, ncol = numFieldCol)
    tempLayoutPlotNum <- matrix(0, nrow = numFieldRow, ncol = numFieldCol)
    tempLayoutRep <- matrix(0, nrow = numFieldRow, ncol = numFieldCol)
    tempLayoutBlock <- matrix(0, nrow = numFieldRow, ncol = numFieldCol)
    
    rowIndexLL1 <- 0
    rowIndexUL1 <- 0
    colIndexLL1 <- 0
    colIndexUL1 <- 0
    
    for (j in 1:numRep){
      
      tempRepTrmt <- matrix(0, nrow = numRowPerRep, ncol = numColPerRep)
      tempRepPlotNum <- matrix(0, nrow = numRowPerRep, ncol = numColPerRep)
      tempRepRep <- matrix(0, nrow = numRowPerRep, ncol = numColPerRep)
      tempRepBlock <- matrix(0, nrow = numRowPerRep, ncol = numColPerRep)
      
      rowIndexLL <- 0
      rowIndexUL <- 0
      colIndexLL <- 0
      colIndexUL <- 0
      
      for(k in 1:numBlk){
        
        tempBlkTrmt <- matrix(0, nrow = numRowPerBlk, ncol = numColPerBlk)
        tempBlkPlotNum <- matrix(0, nrow = numRowPerBlk, ncol = numColPerBlk)
        tempBlkRep <- matrix(0, nrow = numRowPerBlk, ncol = numColPerBlk)
        tempBlkBlock <- matrix(0, nrow = numRowPerBlk, ncol = numColPerBlk)
        
        if(topToBottom & serpentine & numColPerBlk == 1){
          range <- getRange (index = k, topToBottom = topToBottom, serpentine = FALSE,
                             numRowPerGrp = numRowPerBlk, numColPerGrp = numColPerBlk,
                             numGrpRow = numBlkRow, numGrpCol = numBlkCol,
                             numRow = numRowPerRep, numCol = numColPerRep,
                             rowIndexLL = rowIndexLL, rowIndexUL = rowIndexUL,
                             colIndexLL = colIndexLL, colIndexUL = colIndexUL)
        } else if(!topToBottom & serpentine & numRowPerBlk == 1){
          range <- getRange (index = k, topToBottom = topToBottom, serpentine = FALSE,
                             numRowPerGrp = numRowPerBlk, numColPerGrp = numColPerBlk,
                             numGrpRow = numBlkRow, numGrpCol = numBlkCol,
                             numRow = numRowPerRep, numCol = numColPerRep,
                             rowIndexLL = rowIndexLL, rowIndexUL = rowIndexUL,
                             colIndexLL = colIndexLL, colIndexUL = colIndexUL)
        } else{
          range <- getRange (index = k, topToBottom = topToBottom, serpentine = serpentine,
                             numRowPerGrp = numRowPerBlk, numColPerGrp = numColPerBlk,
                             numGrpRow = numBlkRow, numGrpCol = numBlkCol,
                             numRow = numRowPerRep, numCol = numColPerRep,
                             rowIndexLL = rowIndexLL, rowIndexUL = rowIndexUL,
                             colIndexLL = colIndexLL, colIndexUL = colIndexUL)
        }
        
        rowIndexLL <- range$rowIndexLL
        rowIndexUL <- range$rowIndexUL
        colIndexLL <- range$colIndexLL
        colIndexUL <- range$colIndexUL
        
        if(k == 1){
          m <- ((j-1)*numTrmt) + 1
        } else{
          m <- m + blksize
        }
        
        if(topToBottom){
          
          tempBlkTrmt <- matrix(fieldbook[m:(m+blksize-1),"ID"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = FALSE)
          tempBlkPlotNum <- matrix(fieldbook[m:(m+blksize-1),"PlotNumber"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = FALSE)
          tempBlkRep <- matrix(fieldbook[m:(m+blksize-1),"REP"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = FALSE)
          tempBlkBlock <- matrix(fieldbook[m:(m+blksize-1),"Block"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = FALSE)
          
          if (serpentine & numColPerBlk != 1) { 
            for (z in seq(2, numColPerBlk, by = 2)) { 
              tempBlkTrmt[,z] <- rev(tempBlkTrmt[,z]) 
              tempBlkPlotNum[,z] <- rev(tempBlkPlotNum[,z])
              tempBlkRep[,z] <- rev(tempBlkRep[,z])
              tempBlkBlock[,z] <- rev(tempBlkBlock[,z]) }}
          
          tempRepTrmt[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkTrmt
          tempRepPlotNum[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkPlotNum
          tempRepRep[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkRep
          tempRepBlock[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkBlock
          
        } else{
          
          tempBlkTrmt <- matrix(fieldbook[m:(m+blksize-1),"ID"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = TRUE)
          tempBlkPlotNum <- matrix(fieldbook[m:(m+blksize-1),"PlotNumber"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = TRUE)
          tempBlkRep <- matrix(fieldbook[m:(m+blksize-1),"REP"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = TRUE)
          tempBlkBlock <- matrix(fieldbook[m:(m+blksize-1),"Block"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = TRUE)
          
          if (serpentine & numRowPerBlk != 1) { 
            for (z in seq(2, numRowPerBlk, by = 2)) { 
              tempBlkTrmt[z,] <- rev(tempBlkTrmt[z,]) 
              tempBlkPlotNum[z,] <- rev(tempBlkPlotNum[z,])
              tempBlkRep[z,] <- rev(tempBlkRep[z,])
              tempBlkBlock[z,] <- rev(tempBlkBlock[z,]) }}
          
          tempRepTrmt[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkTrmt
          tempRepPlotNum[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkPlotNum
          tempRepRep[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkRep
          tempRepBlock[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkBlock
          
        }
        
      }
      
      if(topToBottom & serpentine & numColPerBlk == 1 & numColPerRep != 1 & numRowPerRep != numFieldRow){
        for (m in seq(2, numColPerRep, by = 2)) { 
          tempRepTrmt[,m] <- rev(tempRepTrmt[,m]) 
          tempRepPlotNum[,m] <- rev(tempRepPlotNum[,m])
          tempRepRep[,m] <- rev(tempRepRep[,m])
          tempRepBlock[,m] <- rev(tempRepBlock[,m]) }
      } else if(!topToBottom & serpentine & numRowPerBlk == 1 & numRowPerRep != 1 & numColPerRep != numFieldCol){
        for (m in seq(2, numRowPerRep, by = 2)) { 
          tempRepTrmt[m,] <- rev(tempRepTrmt[m,]) 
          tempRepPlotNum[m,] <- rev(tempRepPlotNum[m,])
          tempRepRep[m,] <- rev(tempRepRep[m,])
          tempRepBlock[m,] <- rev(tempRepBlock[m,]) }
      }
      
      range <- getRange (index = j, topToBottom = topToBottom, serpentine = serpentine,
                         numRowPerGrp = numRowPerRep, numColPerGrp = numColPerRep,
                         numGrpRow = numRepRow, numGrpCol = numRepCol,
                         numRow = numFieldRow, numCol = numFieldCol,
                         rowIndexLL = rowIndexLL1, rowIndexUL = rowIndexUL1,
                         colIndexLL = colIndexLL1, colIndexUL = colIndexUL1)
      
      rowIndexLL1 <- range$rowIndexLL
      rowIndexUL1 <- range$rowIndexUL
      colIndexLL1 <- range$colIndexLL
      colIndexUL1 <- range$colIndexUL
      
      tempLayoutTrmt[rowIndexLL1:rowIndexUL1,colIndexLL1:colIndexUL1] <- tempRepTrmt
      tempLayoutPlotNum[rowIndexLL1:rowIndexUL1,colIndexLL1:colIndexUL1] <- tempRepPlotNum
      tempLayoutRep[rowIndexLL1:rowIndexUL1,colIndexLL1:colIndexUL1] <- tempRepRep
      tempLayoutBlock[rowIndexLL1:rowIndexUL1,colIndexLL1:colIndexUL1] <- tempRepBlock
    }
    
    if(topToBottom & serpentine & numColPerBlk == 1 & numColPerRep != 1 & numRowPerRep == numFieldRow){
      for (m in seq(2, numFieldCol, by = 2)) { 
        tempLayoutTrmt[,m] <- rev(tempLayoutTrmt[,m]) 
        tempLayoutPlotNum[,m] <- rev(tempLayoutPlotNum[,m])
        tempLayoutRep[,m] <- rev(tempLayoutRep[,m])
        tempLayoutBlock[,m] <- rev(tempLayoutBlock[,m]) }
    } else if(!topToBottom & serpentine & numRowPerBlk == 1 & numRowPerRep != 1 & numColPerRep == numFieldCol){
      for (m in seq(2, numFieldRow, by = 2)) { 
        tempLayoutTrmt[m,] <- rev(tempLayoutTrmt[m,]) 
        tempLayoutPlotNum[m,] <- rev(tempLayoutPlotNum[m,])
        tempLayoutRep[m,] <- rev(tempLayoutRep[m,])
        tempLayoutBlock[m,] <- rev(tempLayoutBlock[m,]) }
    } else if(topToBottom & serpentine & numColPerBlk == 1 & numColPerRep == 1 & numFieldCol != 1){
      for (m in seq(2, numFieldCol, by = 2)) { 
        tempLayoutTrmt[,m] <- rev(tempLayoutTrmt[,m]) 
        tempLayoutPlotNum[,m] <- rev(tempLayoutPlotNum[,m])
        tempLayoutRep[,m] <- rev(tempLayoutRep[,m])
        tempLayoutBlock[,m] <- rev(tempLayoutBlock[,m]) }
    } else if(!topToBottom & serpentine & numRowPerBlk == 1 & numRowPerRep == 1 & numFieldRow != 1){
      for (m in seq(2, numFieldRow, by = 2)) { 
        tempLayoutTrmt[m,] <- rev(tempLayoutTrmt[m,]) 
        tempLayoutPlotNum[m,] <- rev(tempLayoutPlotNum[m,])
        tempLayoutRep[m,] <- rev(tempLayoutRep[m,])
        tempLayoutBlock[m,] <- rev(tempLayoutBlock[m,]) }
    }
    
    tempfbook <- as.data.frame.table(tempLayoutPlotNum)
    tempfbook$Var1 <- as.numeric(tempfbook$Var1)
    tempfbook$Var2 <- as.numeric(tempfbook$Var2)
    
    fieldbook <- merge(fieldbook, tempfbook, by.x = "PlotNumber", by.y = "Freq", all.x = TRUE)
    names(fieldbook)[6:7] <- c("ROW", "RANGE")
    
    trmtLayout <- list(tempLayoutTrmt, tempLayoutPlotNum, tempLayoutRep, tempLayoutBlock)
    
    fieldbook <- fieldbook[,c("ID", "ENTRY", "ROW", "RANGE", "REP", "PlotNumber", "Block")]
  }
  
  return(invisible(list(fieldbook = fieldbook, plan = trmtLayout, designparam = designparam)))
  
}
