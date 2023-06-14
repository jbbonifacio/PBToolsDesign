# -------------------------------------------------------------------------------------
# Name             : randomizeAlphaLatticeWithDiagCheck
# Description      : Generate randomization for Alpha Lattice design (from agricolae)
#                    with diagonal checks
# R Version        : 4.1.2
# -------------------------------------------------------------------------------------
# Author           : Justine B. Dayrit 
# Author Email     : j.bonifacio@irri.org
# Script Version   : 1.2.0
# Modified by      : Justine B. Dayrit
# Date Modified    : 12.01.2022
# Changes made     : updated arrangement of layouts
# -------------------------------------------------------------------------------------
#' @name randomizeAlphaLatticeWithDiagCheck
#' @aliases randomizeAlphaLatticeWithDiagCheck
#  @aliases randomizeAlphaLatticeWithDiagCheck.default
#' @title Randomization for Alpha Lattice design with diagonal checks
#'
#' @description Generate randomization.
#'
#' @param numTrmt number of treatment
#' @param numSpatCheck number of spatial checks
#' @param numRep number of replicates
#' @param numBlk number of blocks per replicate
#' @param blksize number of plots per block
#' @param numRowPerBlk number of row per block
#' @param numColPerBlk number of column per block
#' @param numRowPerRep number of row per replicate
#' @param numColPerRep number of column per replicate
#' @param numFieldRow number of field rows
#' @param numFieldCol number of field columns
#' @param topToBottom a logical variable indicating whether plot number will be written from to to bottom
#' @param serpentine a logical variable indicating whether plot number will be arranged as serpentine order
#' @param trmtList NULL or vector or character containing the levels of the treatment
#' 
#' @return A list containing dataframe and statistical design array.
#'
#' @examples
#' randomizeAlphaLatticeWithDiagCheck(numTrmt = 27, numSpatCheck = 1, 
#'                                    numFieldRow = 9, numFieldCol = 6, 
#'                                    numRep = 2, numRowPerRep = 9, numColPerRep = 3, 
#'                                    numBlk = 9, blksize = 3, numRowPerBlk = 3, numColPerBlk = 1,
#'                                    topToBottom = TRUE, serpentine = TRUE, 
#'                                    trmtList = c(paste("Entry",1:27,sep=""),paste("SpatialCheck",1,sep="")))
#' 
# -------------------------------------------------------------------------------------

# Parameters for the alpha design: 
# I. r=2, k <= s; 
# II.	r=3, s odd, k <= s; 
# III.r=3, s even, k <= s-1; 
# IV.	r=4, s odd but not a multiple of 3, k<=s

# r = replications
# s = number of blocks 
# k = size of block 
# Number of treatment is equal to k*s

randomizeAlphaLatticeWithDiagCheck <- function(numTrmt, numSpatCheck, numFieldRow, numFieldCol,
                                               numRep, numRowPerRep, numColPerRep,
                                               numBlk, blksize, numRowPerBlk, numColPerBlk,
                                               topToBottom = NULL, serpentine = NULL, trmtList = NULL) {
  
  sink(tempfile())
  result <- try(temp <- design.alpha(trt = 1:numTrmt,
                                     k = blksize,
                                     r = numRep,
                                     serie = nchar(numTrmt)),
                silent = TRUE)
  sink()
  
  if (all(class(result) == "try-error")) {
    stop(paste("Error found in using design.alpha from agricolae package"))
  }
  
  if(numRowPerBlk == 1){
    numBlkPerRow <- numFieldCol/numColPerBlk
    numFieldCol <- numFieldCol + numBlkPerRow
    numColPerRep <- numColPerRep + (numColPerRep/numColPerBlk)
    numColPerBlk <- numColPerBlk + 1
    percentCheck <- (numBlkPerRow*numFieldRow)/(numFieldCol*numFieldRow)
    
    dlayout <- diag.layout(numFieldCol, numFieldRow, percentCheck, plot = FALSE)
    
    temptrmt <- as.vector(t(dlayout))
  } else if(numRowPerBlk == blksize){
    numBlkPerCol <- numFieldRow/numRowPerBlk
    numFieldRow <- numFieldRow + numBlkPerCol
    numRowPerRep <- numRowPerRep + (numRowPerRep/numRowPerBlk)
    numRowPerBlk <- numRowPerBlk + 1
    percentCheck <- (numBlkPerCol*numFieldCol)/(numFieldRow*numFieldCol)
    
    dlayout <- diag.layout(numFieldRow, numFieldCol, percentCheck, plot = FALSE)
    dlayout <- t(dlayout)
    
    temptrmt <- as.vector(dlayout)
  } else {
    stop(paste("This feature is not yet available"))
  }
  
  numRepCol <- numFieldCol/numColPerRep 
  numRepRow <- numFieldRow/numRowPerRep
  numBlkRow <- numRowPerRep/numRowPerBlk 
  numBlkCol <- numColPerRep/numColPerBlk
  blksize <- blksize + 1
  
  spatCheckNumS <- numTrmt+1
  spatCheckNumE <- numTrmt+numSpatCheck
  
  if(sum(dlayout)%%numSpatCheck == 0){
    tempSpat <- rep(spatCheckNumS:spatCheckNumE,each = (sum(dlayout)/numSpatCheck))
    tempSpatOrder <- data.frame(tempSpat , tempPlotNum = sample(length(tempSpat), length(tempSpat), replace = FALSE))
    tempSpatOrder <- tempSpatOrder[order(tempSpatOrder[,"tempPlotNum"]),]
  } else{
    tempSpat <- c(rep(spatCheckNumS:spatCheckNumE,each = floor(sum(dlayout)/numSpatCheck)),sample(spatCheckNumS:spatCheckNumE, (sum(dlayout)%%numSpatCheck), replace = FALSE))
    tempSpatOrder <- data.frame(tempSpat , tempPlotNum = sample(length(tempSpat), length(tempSpat), replace = FALSE))
    tempSpatOrder <- tempSpatOrder[order(tempSpatOrder[,"tempPlotNum"]),]
  }
  
  temptrmt[temptrmt == 1] <- as.numeric(tempSpatOrder$tempSpat)
  temptrmt[temptrmt == 0] <- as.numeric(as.vector(result$book$`1:numTrmt`))

  fieldbook <- data.frame("PlotNumber" = 1:length(dlayout), "Block" = rep(rep(1:numBlk,each = blksize),numRep), "ENTRY" = temptrmt, "REP" = rep(1:numRep, each = length(dlayout)/numRep))
  fieldbook$PlotNumber <- as.numeric(paste(fieldbook$REP, paste(rep(0,nchar(length(dlayout)/numRep)), collapse=""), sep = "")) + 1:(length(dlayout)/numRep)
  
  #merge treatment list with fielbook by entry number
  trmtList <- data.frame("ENTRY" = 1:(numTrmt+numSpatCheck), "ID" = trmtList)
  fieldbook <- merge(fieldbook, trmtList, all.x = TRUE)
  fieldbook <- fieldbook[order(as.numeric(fieldbook$PlotNumber)),]
  rownames(fieldbook) <- 1:nrow(fieldbook)
    
  fieldbook$REP <- as.numeric(fieldbook$REP)
  
  if(topToBottom){
    if(numRowPerBlk == blksize){
      tempLayoutTrmt <- matrix(fieldbook$ID, nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)  
    } else{
      tempLayoutTrmt <- matrix(fieldbook$ID, nrow = numFieldRow, ncol = numFieldCol, byrow = TRUE)  
    }
  } else{
    if(numRowPerBlk == 1){
      tempLayoutTrmt <- matrix(fieldbook$ID, nrow = numFieldRow, ncol = numFieldCol, byrow = TRUE) 
    } else{
      tempLayoutTrmt <- matrix(fieldbook$ID, nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)  
    }
  }
  
  tempLayoutPlotNum <- matrix(0, nrow = numFieldRow, ncol = numFieldCol)
  tempLayoutRep <- matrix(0, nrow = numFieldRow, ncol = numFieldCol)
  tempLayoutBlock <- matrix(0, nrow = numFieldRow, ncol = numFieldCol)
  
  rowIndexLL1 <- 0
  rowIndexUL1 <- 0
  colIndexLL1 <- 0
  colIndexUL1 <- 0
  
  for (j in 1:numRep){
    
    tempRepPlotNum <- matrix(0, nrow = numRowPerRep, ncol = numColPerRep)
    tempRepRep <- matrix(0, nrow = numRowPerRep, ncol = numColPerRep)
    tempRepBlock <- matrix(0, nrow = numRowPerRep, ncol = numColPerRep)
    
    rowIndexLL <- 0
    rowIndexUL <- 0
    colIndexLL <- 0
    colIndexUL <- 0
    
    for(k in 1:numBlk){
      
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
        m <- ((j-1)*(length(dlayout)/numRep)) + 1
      } else{
        m <- m + blksize
      }
      
      if(topToBottom){
        
        tempBlkPlotNum <- matrix(fieldbook[m:(m+blksize-1),"PlotNumber"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = FALSE)
        tempBlkRep <- matrix(fieldbook[m:(m+blksize-1),"REP"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = FALSE)
        tempBlkBlock <- matrix(fieldbook[m:(m+blksize-1),"Block"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = FALSE)
        
        if (serpentine & numColPerBlk != 1) { 
          for (z in seq(2, numColPerBlk, by = 2)) { 
            tempBlkPlotNum[,z] <- rev(tempBlkPlotNum[,z])
            tempBlkRep[,z] <- rev(tempBlkRep[,z])
            tempBlkBlock[,z] <- rev(tempBlkBlock[,z]) }}
        
        tempRepPlotNum[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkPlotNum
        tempRepRep[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkRep
        tempRepBlock[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkBlock
        
      } else{
        
        tempBlkPlotNum <- matrix(fieldbook[m:(m+blksize-1),"PlotNumber"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = TRUE)
        tempBlkRep <- matrix(fieldbook[m:(m+blksize-1),"REP"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = TRUE)
        tempBlkBlock <- matrix(fieldbook[m:(m+blksize-1),"Block"], nrow = numRowPerBlk, ncol = numColPerBlk, byrow = TRUE)
        
        if (serpentine & numRowPerBlk != 1) { 
          for (z in seq(2, numRowPerBlk, by = 2)) { 
            tempBlkPlotNum[z,] <- rev(tempBlkPlotNum[z,])
            tempBlkRep[z,] <- rev(tempBlkRep[z,])
            tempBlkBlock[z,] <- rev(tempBlkBlock[z,]) }}
        
        tempRepPlotNum[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkPlotNum
        tempRepRep[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkRep
        tempRepBlock[rowIndexLL:rowIndexUL,colIndexLL:colIndexUL] <- tempBlkBlock
        
      }
      
    }
    
    if(topToBottom & serpentine & numColPerBlk == 1 & numColPerRep != 1 & numRowPerRep != numFieldRow){
      for (m in seq(2, numColPerRep, by = 2)) { 
        tempRepPlotNum[,m] <- rev(tempRepPlotNum[,m])
        tempRepRep[,m] <- rev(tempRepRep[,m])
        tempRepBlock[,m] <- rev(tempRepBlock[,m]) }
    } else if(!topToBottom & serpentine & numRowPerBlk == 1 & numRowPerRep != 1 & numColPerRep != numFieldCol){
      for (m in seq(2, numRowPerRep, by = 2)) { 
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
    
    tempLayoutPlotNum[rowIndexLL1:rowIndexUL1,colIndexLL1:colIndexUL1] <- tempRepPlotNum
    tempLayoutRep[rowIndexLL1:rowIndexUL1,colIndexLL1:colIndexUL1] <- tempRepRep
    tempLayoutBlock[rowIndexLL1:rowIndexUL1,colIndexLL1:colIndexUL1] <- tempRepBlock
  
  }
  
  # if(topToBottom & serpentine & numColPerBlk == 1 & numColPerRep != 1 & numRowPerRep == numFieldRow){
  #   for (m in seq(2, numFieldCol, by = 2)) { 
  #     tempLayoutPlotNum[,m] <- rev(tempLayoutPlotNum[,m])
  #     tempLayoutRep[,m] <- rev(tempLayoutRep[,m])
  #     tempLayoutBlock[,m] <- rev(tempLayoutBlock[,m]) }
  # } else if(!topToBottom & serpentine & numRowPerBlk == 1 & numRowPerRep != 1 & numColPerRep == numFieldCol){
  #   for (m in seq(2, numFieldRow, by = 2)) { 
  #     tempLayoutPlotNum[m,] <- rev(tempLayoutPlotNum[m,])
  #     tempLayoutRep[m,] <- rev(tempLayoutRep[m,])
  #     tempLayoutBlock[m,] <- rev(tempLayoutBlock[m,]) }
  # } else if(topToBottom & serpentine & numColPerBlk == 1 & numColPerRep == 1 & numFieldCol != 1){
  #   for (m in seq(2, numFieldCol, by = 2)) { 
  #     tempLayoutPlotNum[,m] <- rev(tempLayoutPlotNum[,m])
  #     tempLayoutRep[,m] <- rev(tempLayoutRep[,m])
  #     tempLayoutBlock[,m] <- rev(tempLayoutBlock[,m]) }
  # } else if(!topToBottom & serpentine & numRowPerBlk == 1 & numRowPerRep == 1 & numFieldRow != 1){
  #   for (m in seq(2, numFieldRow, by = 2)) { 
  #     tempLayoutPlotNum[m,] <- rev(tempLayoutPlotNum[m,])
  #     tempLayoutRep[m,] <- rev(tempLayoutRep[m,])
  #     tempLayoutBlock[m,] <- rev(tempLayoutBlock[m,]) }
  # }
  
  tempfbook <- as.data.frame.table(tempLayoutPlotNum)
  tempfbook$Var1 <- as.numeric(tempfbook$Var1)
  tempfbook$Var2 <- as.numeric(tempfbook$Var2)
  
  tempfbook <- tempfbook[order(as.numeric(tempfbook$Freq)),]
  rownames(tempfbook) <- 1:nrow(tempfbook)
  
  if(topToBottom){
    if(numRowPerBlk == blksize){
      fieldbook$PlotNumber <- as.vector(tempLayoutPlotNum)
      fieldbook$REP <- as.vector(tempLayoutRep)
      fieldbook$Block <- as.vector(tempLayoutBlock) 
    } else{
      fieldbook$PlotNumber <- as.vector(t(tempLayoutPlotNum))
      fieldbook$REP <- as.vector(t(tempLayoutRep))
      fieldbook$Block <- as.vector(t(tempLayoutBlock))
    }
  } else{
    if(numRowPerBlk == 1){
      fieldbook$PlotNumber <- as.vector(t(tempLayoutPlotNum))
      fieldbook$REP <- as.vector(t(tempLayoutRep))
      fieldbook$Block <- as.vector(t(tempLayoutBlock))
    } else{
      fieldbook$PlotNumber <- as.vector(tempLayoutPlotNum)
      fieldbook$REP <- as.vector(tempLayoutRep)
      fieldbook$Block <- as.vector(tempLayoutBlock) 
    }
  }
  
  fieldbook <- fieldbook[order(as.numeric(fieldbook$PlotNumber)),]
  rownames(fieldbook) <- 1:nrow(fieldbook)
  
  fieldbook <- merge(fieldbook, tempfbook, by.x = "PlotNumber", by.y = "Freq", all.x = TRUE)
  names(fieldbook)[6:7] <- c("ROW", "RANGE")
  
  trmtLayout <- list(tempLayoutTrmt, tempLayoutPlotNum, tempLayoutRep, tempLayoutBlock)
  
  fieldbook <- fieldbook[,c("ID", "ENTRY", "ROW", "RANGE", "REP", "PlotNumber", "Block")]
  
  return(invisible(list(fieldbook = fieldbook, plan = trmtLayout)))
  
}
