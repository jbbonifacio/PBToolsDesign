# -------------------------------------------------------------------------------------
# Name             : layoutAlphaLatticeWithRepCheck
# Description      : Generate layout for Alpha Lattice Design with Repeated Checks
# R Version        : 4.0.3
# -------------------------------------------------------------------------------------
# Author           : Justine B. Bonifacio
# Author Email     : j.bonifacio@irri.org
# Script Version   : 1.0.0
# Modified by      : 
# Date Modified    : 
# Changes made     : 
# -------------------------------------------------------------------------------------
#' @name layoutAlphaLatticeWithRepCheck
#' @aliases layoutAlphaLatticeWithRepCheck
#  @aliases layoutAlphaLatticeWithRepCheck.default
#' @title Layout Creation for Alpha Lattice Design with Repeated Checks
#' 
#' @description Creates layout for Alpha Lattice Design with Repeated Checks.
#'
#' @param fieldbook a dataframe which is a result from the function \code{randomizeAlphaLatticeWithRepCheck} 
#' @param trmtLayout a list which is a result from the function \code{randomizeAlphaLatticeWithRepCheck} 
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

layoutAlphaLatticeWithRepCheck <- function(fieldbook, trmtLayout, numFieldRow, numFieldCol,
                                           numRowPerRep, numColPerRep, numBlk,
                                           numRowPerBlk, numColPerBlk,
                                           serpentine = FALSE, topToBottom = TRUE) { #UseMethod("layoutAlphaLatticeWithRepCheck") }
  
#layoutAlphaLatticeWithRepCheck.default <- function(fieldbook, trmtLayout,
#                                                   numFieldRow, numFieldCol,
#                                                   numRowPerRep, numColPerRep,
#                                                   numBlk, numRowPerBlk, numColPerBlk,
#                                                   serpentine = FALSE,
#                                                   topToBottom = TRUE) {
  
  numRepRow <- numFieldRow/numRowPerRep
  numRepCol <- numFieldCol/numColPerRep
  
  fieldbook$Trial <- as.factor(fieldbook$Trial)
  fieldbook$ID <- as.factor(fieldbook$ID)
  fieldbook$REP <- as.factor(fieldbook$REP)
  
  numBlkRow <- numRowPerRep/numRowPerBlk
  numBlkCol <- numColPerRep/numColPerBlk
  blksize <- numRowPerBlk*numColPerBlk
  
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
    plotNum <- matrix(as.numeric(paste(fieldbook$REP, paste(rep(0,nchar(numRowPerRep*numColPerRep)), collapse = ""), sep = "")), 
                      nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)
    plotNumPerBlk <- matrix(1:blksize, nrow = numRowPerBlk, ncol = numColPerBlk, byrow = FALSE)
    if (serpentine & numColPerBlk != 1) { for (z in seq(2, numColPerBlk, by = 2)) { plotNumPerBlk[,z] <- rev(plotNumPerBlk[,z]) }} # serpentine per block
  } else{
    fieldbook$REP <- as.vector(repNum)
    plotNum <- matrix(as.numeric(paste(fieldbook$REP, paste(rep(0,nchar(numRowPerRep*numColPerRep)), collapse = ""), sep = "")), 
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
    
    range <- getRange(index = j, topToBottom = topToBottom, serpentine = serpentine,
                      numRowPerGrp = numRowPerBlk, numColPerGrp = numColPerBlk,
                      numGrpRow = numBlkRow, numGrpCol = numBlkCol,
                      numRow = numRowPerRep, numCol = numColPerRep,
                      rowIndexLL = rowIndexLL, rowIndexUL = rowIndexUL,
                      colIndexLL = colIndexLL, colIndexUL = colIndexUL)
    
    rowIndexLL <- range$rowIndexLL
    rowIndexUL <- range$rowIndexUL
    colIndexLL <- range$colIndexLL
    colIndexUL <- range$colIndexUL
    
    blkNumPerRep[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tempBlkNumPerRep
    plotNumPerRep[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tempPlotNumPerRep
    tempPlotNumPerRep <- tempPlotNumPerRep + blksize
    
  }
  
  plotNum <- writePlotNumInfo(info = plotNumPerRep, table =  plotNum, tableGrp = NULL, 
                              numGrpRow = numRepRow, numGrpCol = numRepCol, 
                              numRowPerGrp = numRowPerRep, numColPerGrp = numColPerRep, incrementWithin = FALSE,
                              topToBottom = topToBottom)[[1]]
  
  # block information
  blkNum <- writePlotNumInfo(info = blkNumPerRep, table =  blkNum, tableGrp = NULL, 
                             numGrpRow = numRepRow, numGrpCol = numRepCol, 
                             numRowPerGrp = numRowPerRep, numColPerGrp = numColPerRep, incrementWithin = FALSE,
                             topToBottom = topToBottom)[[1]]
  
  grpInfo <- merge(as.data.frame.table(repNum), as.data.frame.table(blkNum), 
                   by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
  grpInfo <- merge(grpInfo, as.data.frame.table(plotNum),by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
  names(grpInfo)[3:ncol(grpInfo)] <- c("Rep", "Block", "PlotNumber")  
  
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
  
  newfbook <- merge(fieldbook, tmpfbook, by.x = c("Trial", "ID", "ROW", "RANGE"), by.y = c("Trial", "Freq", "Var1", "Var2"))  
  exptLayout[[4]] <- blkNum
  names(exptLayout) <- c("TrmtLayout", "PlotNumLayout", "RepLayout", "BlockLayout")
  
  return(list(fieldbook = newfbook, plan = exptLayout))
  
} ## end of layoutAlphaLatticeWithRepCheck function
