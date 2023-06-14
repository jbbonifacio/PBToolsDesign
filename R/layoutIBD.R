# -------------------------------------------------------------------------------------
# Name             : layoutIBD 
# Description      : Generate layout for Incomplete Block Designs
# R Version        : 4.0.1 
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 1.3.1
# Modified by      : Justine B. Dayrit
# Date Modified    : 2022.10.07
# Changes made     : updated matrix data to include trial 1
# -------------------------------------------------------------------------------------
#' @name layoutIBD
#' @aliases layoutIBD
#  @aliases layoutIBD.default
#' @title Layout Creation for Incomplete Block Designs
#'
#' @description Creates layout Incomplete Block Designs including Alpha Lattice and Row-Column Design.
#'
#' @param fieldbook a dataframe which is a result from the function \code{randomizeIBD} 
#' @param trmtLayout a list which is a result from the function \code{randomizeIBD} 
#' @param numFieldRow number of field rows
#' @param numFieldCol number of field columns
#' @param numRowPerRep number of rows per replicate
#' @param numColPerRep number of columns per replicate
#' @param numBlk number of blocks
#' @param numRowPerBlk number of rows per block
#' @param numColPerBlk number of columns per block
#' @param topToBottom logical, whether plot number will be written from to to bottom
#' @param serpentine logical, whether plot number will be arranged as serpentine order
#' 
#' @return list containing a fieldbook, the layout and plot numbers.
#' 
# -------------------------------------------------------------------------------------

layoutIBD <- function(fieldbook, trmtLayout,
                      numFieldRow, numFieldCol,
                      numRowPerRep, numColPerRep, numBlk = NULL,
                      numRowPerBlk = NULL, numColPerBlk = NULL,
                      serpentine = FALSE, topToBottom = TRUE) { #UseMethod("layoutIBD") }
  
#layoutIBD.default <- function(fieldbook, trmtLayout, 
#                              numFieldRow, numFieldCol,
#                              numRowPerRep, numColPerRep, 
#                              numBlk = NULL, numRowPerBlk = NULL, numColPerBlk = NULL,
#                              serpentine = FALSE,
#                              topToBottom = TRUE) {
  
  numRepRow <- numFieldRow/numRowPerRep
  numRepCol <- numFieldCol/numColPerRep

  if (!is.null(numRowPerBlk) && !is.null(numColPerBlk)) {
    # Design == Alpha Lattice

    fieldbook$Trial <- as.factor(fieldbook$Trial)
    fieldbook$ID <- as.factor(fieldbook$ID)
    fieldbook$REP <- as.factor(fieldbook$REP)

    numBlkRow <- numRowPerRep/numRowPerBlk
    numBlkCol <- numColPerRep/numColPerBlk
    blksize <- length(unique(fieldbook$ID))/numBlk
    
    plotNumPerRep <- matrix(0, nrow = numRowPerRep, ncol = numColPerRep)
    blkNum <- matrix(0, nrow = numFieldRow, ncol = numFieldCol)
    blkNumPerRep <- matrix(0, nrow = numRowPerRep, ncol = numColPerRep)
    
    rowIndexLL <- 0
    rowIndexUL <- 0
    colIndexLL <- 0
    colIndexUL <- 0
    
    #assignment of rep number
    repNum <- matrix(0, nrow = numFieldRow, ncol = numFieldCol)
    for (j in 1:nlevels(fieldbook$REP)){
      tempRep <- matrix(j, nrow = numRowPerRep, ncol = numColPerRep)
      
      range <- getRange(index = j, topToBottom = topToBottom, serpentine = serpentine,
                        numRowPerGrp = numRowPerRep, numColPerGrp = numColPerRep,
                        numGrpRow = numRepRow, numGrpCol = numRepCol,
                        numRow = numFieldRow, numCol = numFieldCol,
                        rowIndexLL = rowIndexLL, rowIndexUL = rowIndexUL,
                        colIndexLL = colIndexLL, colIndexUL = colIndexUL)

      rowIndexLL <- range$rowIndexLL
      rowIndexUL <- range$rowIndexUL
      colIndexLL <- range$colIndexLL
      colIndexUL <- range$colIndexUL
      
      repNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tempRep

    }
    
    if (topToBottom){
      fieldbook$REP <- as.vector(repNum)
      plotNum <- matrix(as.numeric(paste(fieldbook[fieldbook$Trial == 1,"REP"], paste(rep(0,nchar(nlevels(fieldbook$ID))), collapse = ""), sep = "")), 
                        nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)
      plotNumPerBlk <- matrix(1:blksize, nrow = numRowPerBlk, ncol = numColPerBlk, byrow = FALSE)
      if (serpentine & numColPerBlk != 1) { for (z in seq(2, numColPerBlk, by = 2)) { plotNumPerBlk[,z] <- rev(plotNumPerBlk[,z]) }} # serpentine per block
    } else{
      fieldbook$REP <- as.vector(repNum)
      plotNum <- matrix(as.numeric(paste(fieldbook[fieldbook$Trial == 1,"REP"], paste(rep(0,nchar(nlevels(fieldbook$ID))), collapse = ""), sep = "")), 
                        nrow = numFieldRow, ncol = numFieldCol)
      plotNumPerBlk <- matrix(1:blksize, nrow = numRowPerBlk, ncol = numColPerBlk, byrow = TRUE)
      if (serpentine & numRowPerBlk != 1) { for (z in seq(2, numRowPerBlk, by = 2)) { plotNumPerBlk[z,] <- rev(plotNumPerBlk[z,]) }} # serpentine per block
    }
    
    #assignment of block and plot number
    blkNumPerRep <- matrix(0, nrow = numRowPerRep, ncol = numColPerRep)
    plotNumPerRep <- matrix(0, nrow = numRowPerRep, ncol = numColPerRep)
    tempPlotNumPerRep <- plotNumPerBlk
    
    for (j in 1:numBlk){
      tempBlkNumPerRep <- matrix(j, nrow = numRowPerBlk, ncol = numColPerBlk)
    
      if(topToBottom & serpentine & numColPerBlk == 1){
        range <- getRange(index = j, topToBottom = topToBottom, serpentine = FALSE,
                          numRowPerGrp = numRowPerBlk, numColPerGrp = numColPerBlk,
                          numGrpRow = numBlkRow, numGrpCol = numBlkCol,
                          numRow = numRowPerRep, numCol = numColPerRep,
                          rowIndexLL = rowIndexLL, rowIndexUL = rowIndexUL,
                          colIndexLL = colIndexLL, colIndexUL = colIndexUL)
      } else if (!topToBottom & serpentine & numRowPerBlk == 1){
        range <- getRange(index = j, topToBottom = topToBottom, serpentine = FALSE,
                          numRowPerGrp = numRowPerBlk, numColPerGrp = numColPerBlk,
                          numGrpRow = numBlkRow, numGrpCol = numBlkCol,
                          numRow = numRowPerRep, numCol = numColPerRep,
                          rowIndexLL = rowIndexLL, rowIndexUL = rowIndexUL,
                          colIndexLL = colIndexLL, colIndexUL = colIndexUL)
      } else{
        range <- getRange(index = j, topToBottom = topToBottom, serpentine = serpentine,
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
      
      blkNumPerRep[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tempBlkNumPerRep
      plotNumPerRep[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tempPlotNumPerRep
      tempPlotNumPerRep <- tempPlotNumPerRep + blksize
    
    }
    
    if(topToBottom & serpentine & numColPerBlk == 1 & numColPerRep != 1 & numRowPerRep != numFieldRow){
      for (m in seq(2, numColPerRep, by = 2)) { 
        plotNumPerRep[,m] <- rev(plotNumPerRep[,m]) 
        blkNumPerRep[,m] <- rev(blkNumPerRep[,m]) }
    } else if(!topToBottom & serpentine & numRowPerBlk == 1 & numRowPerRep != 1 & numColPerRep != numFieldCol){
      for (m in seq(2, numRowPerRep, by = 2)) { 
        plotNumPerRep[m,] <- rev(plotNumPerRep[m,]) 
        blkNumPerRep[m,] <- rev(blkNumPerRep[m,]) }
    }
    
    plotNum <- writePlotNumInfo(info = plotNumPerRep, table =  plotNum, tableGrp = NULL, 
                                numGrpRow = numRepRow, numGrpCol = numRepCol, 
                                numRowPerGrp = numRowPerRep, numColPerGrp = numColPerRep, incrementWithin = FALSE,
                                topToBottom = topToBottom)[[1]]
    
    # block information, if design == AlphaLattice
    blkNum <- writePlotNumInfo(info = blkNumPerRep, table =  blkNum, tableGrp = NULL, 
                               numGrpRow = numRepRow, numGrpCol = numRepCol, 
                               numRowPerGrp = numRowPerRep, numColPerGrp = numColPerRep, incrementWithin = FALSE,
                               topToBottom = topToBottom)[[1]]
    
    if(topToBottom & serpentine & numColPerBlk == 1){
      if(numColPerRep != 1){
        if(numRowPerRep == numFieldRow){
          for (m in seq(2, numFieldCol, by = 2)) { 
            plotNum[,m] <- rev(plotNum[,m]) 
            blkNum[,m] <- rev(blkNum[,m]) }
        }
      } else {
        if(numFieldCol != 1){
          for (m in seq(2, numFieldCol, by = 2)) { 
            plotNum[,m] <- rev(plotNum[,m]) 
            blkNum[,m] <- rev(blkNum[,m]) }
        }
      }
    } else if(!topToBottom & serpentine & numRowPerBlk == 1){
      if(numRowPerRep != 1){
        if(numColPerRep == numFieldCol){
          for (m in seq(2, numFieldRow, by = 2)) { 
            plotNum[m,] <- rev(plotNum[m,]) 
            blkNum[m,] <- rev(blkNum[m,]) }
        }
      } else {
        if(numFieldRow != 1){
          for (m in seq(2, numFieldRow, by = 2)) { 
            plotNum[m,] <- rev(plotNum[m,]) 
            blkNum[m,] <- rev(blkNum[m,]) }
        }
      }
    }
  } else {
    # Design == Row-Column
    
    fieldbook$Trial <- as.factor(fieldbook$Trial)
    fieldbook$ID <- as.factor(fieldbook$ID)
    fieldbook$REP <- as.factor(fieldbook$REP)
    fieldbook$RowBlock <- as.factor(fieldbook$RowBlock)
    fieldbook$ColumnBlock <- as.factor(fieldbook$ColumnBlock)
    
    rowIndexLL <- 0
    rowIndexUL <- 0
    colIndexLL <- 0
    colIndexUL <- 0
    
    repNum <- matrix(0, nrow = numFieldRow, ncol = numFieldCol)
    
    for (j in 1:nlevels(fieldbook$REP)){
      tempRep <- matrix(j, nrow = numRowPerRep, ncol = numColPerRep)
      
      range <- getRange(index = j, topToBottom = topToBottom, serpentine = serpentine,
                        numRowPerGrp = numRowPerRep, numColPerGrp = numColPerRep,
                        numGrpRow = numRepRow, numGrpCol = numRepCol,
                        numRow = numFieldRow, numCol = numFieldCol,
                        rowIndexLL = rowIndexLL, rowIndexUL = rowIndexUL,
                        colIndexLL = colIndexLL, colIndexUL = colIndexUL)
      
      rowIndexLL <- range$rowIndexLL
      rowIndexUL <- range$rowIndexUL
      colIndexLL <- range$colIndexLL
      colIndexUL <- range$colIndexUL
      
      repNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tempRep
    }
    
    if (topToBottom) {
      fieldbook$REP <- as.vector(repNum)
      plotNum <- matrix(as.numeric(paste(fieldbook[fieldbook$Trial == 1,"REP"], paste(rep(0,nchar(nlevels(fieldbook$ID))), collapse = ""), sep = "")), 
                        nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)
      plotNumPerRep <- matrix(1:nlevels(fieldbook$ID), nrow = numRowPerRep, ncol = numColPerRep)
      if (serpentine) { for (z in seq(2, numColPerRep, by = 2)) { plotNumPerRep[,z] <- rev(plotNumPerRep[,z]) }} # serpentine per rep
    } else {
      fieldbook$REP <- as.vector(repNum)
      plotNum <- matrix(as.numeric(paste(fieldbook[fieldbook$Trial == 1,"REP"], paste(rep(0,nchar(nlevels(fieldbook$ID))), collapse = ""), sep = "")), 
                        nrow = numFieldRow, ncol = numFieldCol)
      plotNumPerRep <- matrix(1:nlevels(fieldbook$ID), nrow = numRowPerRep, ncol = numColPerRep, byrow = TRUE)
      if (serpentine) { for (z in seq(2, numRowPerRep, by = 2)) { plotNumPerRep[z,] <- rev(plotNumPerRep[z,]) }} # serpentine per rep
    }
    
    rowBlkNum <- matrix(1:numRowPerRep, nrow = numFieldRow, ncol = numFieldCol)
    colBlkNum <- matrix(rep(1:numColPerRep, each = numFieldRow), nrow = numFieldRow, ncol = numFieldCol)

    plotNum <- writePlotNumInfo(info = plotNumPerRep, table =  plotNum, tableGrp = NULL, 
                                # plotNumOrient = plotNumOrient, 
                                numGrpRow = numRepRow, numGrpCol = numRepCol, 
                                numRowPerGrp = numRowPerRep, numColPerGrp = numColPerRep, 
                                topToBottom = topToBottom, incrementWithin = FALSE)[[1]]
    
  } ## end stmt -- if (!is.null(numRowPerBlk) && !is.null(numColPerBlk)) - else stmt

  if (!is.null(numRowPerBlk) && !is.null(numColPerBlk)) {
    grpInfo <- merge(as.data.frame.table(repNum), as.data.frame.table(blkNum), 
                     by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
    grpInfo <- merge(grpInfo, as.data.frame.table(plotNum),by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
    names(grpInfo)[3:ncol(grpInfo)] <- c("Rep", "Block", "PlotNumber")  
  } else {
    grpInfo <- merge(as.data.frame.table(repNum), as.data.frame.table(rowBlkNum), 
                     by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
    grpInfo <- merge(grpInfo, as.data.frame.table(colBlkNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
    grpInfo <- merge(grpInfo, as.data.frame.table(plotNum), 
                     by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"), 
                     suffixes = c(".1", ".2"))
    names(grpInfo)[3:ncol(grpInfo)] <- c("Rep", "RowBlock", "ColumnBlock", "PlotNumber")
  } ## end stmt -- if (!is.null(numRowPerBlk) && !is.null(numColPerBlk)) - else stmt
    
  for (i in (1:length(trmtLayout))) {
    if (i == 1) { tmpfbook <- data.frame(Trial = i, merge(grpInfo, as.data.frame.table(trmtLayout[[i]]), 
                                              by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2")))
    } else {
      tmpfbook <- rbind(tmpfbook,
                        cbind(Trial = i, merge(grpInfo, as.data.frame.table(trmtLayout[[i]]), 
                                               by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))))
    } ## end stmt -- if (i == 1) - else stmt 
  } ## end stmt -- for (i in (1:length(trmtLayout)))
  
  tmpfbook[,"Var1"] <- as.numeric(tmpfbook[,"Var1"])
  tmpfbook[,"Var2"] <- as.numeric(tmpfbook[,"Var2"])
  
  exptLayout <- list()
  exptLayout[[1]] <- trmtLayout
  names(exptLayout[[1]]) <- paste("Trial", 1:length(trmtLayout), sep = "")
  exptLayout[[2]] <- plotNum
  exptLayout[[3]] <- repNum
  
  if (!is.null(numRowPerBlk) && !is.null(numColPerBlk)) {
    newfbook <- merge(fieldbook, tmpfbook, by.x = c("Trial", "ID", "ROW", "RANGE"), by.y = c("Trial", "Freq", "Var1", "Var2"))  
    exptLayout[[4]] <- blkNum
    names(exptLayout) <- c("TrmtLayout", "PlotNumLayout", "RepLayout", "BlockLayout")
  } else {
    newfbook <- merge(fieldbook, tmpfbook, 
                      by.x = c("Trial", "ID", "ROW", "RANGE", "RowBlock", "ColumnBlock"), 
                      by.y = c("Trial", "Freq", "Var1", "Var2","RowBlock", "ColumnBlock"))
    exptLayout[[4]] <- rowBlkNum
    exptLayout[[5]] <- colBlkNum
    names(exptLayout) <- c("TrmtLayout", "PlotNumLayout", "RepLayout", "RowBlockLayout","ColumnBlockLayout")
  } ## end stmt -- if (!is.null(numRowPerBlk) && !is.null(numColPerBlk)) - else stmt
  
  return(list(fieldbook = newfbook, plan = exptLayout))
  
} ## end of layoutIBD function
