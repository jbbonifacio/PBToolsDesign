# -------------------------------------------------------------------------------------
# Name             : randomizeAugmentedAL
# Description      : Generate randomization for Augmented Alpha Lattice designs
# R Version        : 4.1.0
# -------------------------------------------------------------------------------------
# Author           : Justine B. Bonifacio
# Author Email     : j.bonifacio@irri.org
# Script Version   : 1.0.0
# Modified by      : 
# Date Modified    : 
# Changes made     : 
# -------------------------------------------------------------------------------------
#' @name randomizeAugmentedAL
#' @aliases randomizeAugmentedAL
#  @aliases randomizeAugmentedAL.default
#' @title Randomization for Augmented Alpha Lattice designs
#'
#' @description Generate randomization.
#'
#' @param numCheck number of replicated treatments
#' @param numTest number of unreplicated treatments
#' @param numBlk number of blocks per replicate
#' @param numRep number of replicates
#' @param genLayout logical, whether a layout of the design will be generated or not
#' @param numRowPerBlk number of rows per block, if genLayout is TRUE
#' @param numColPerBlk number of columns per block, if genLayout is TRUE
#' @param numRowPerRep number of rows per replicate, if genLayout is TRUE
#' @param numColPerRep number of columns per replicate, if genLayout is TRUE
#' @param numFieldRow number of field rows, if genLayout is TRUE
#' @param numFieldCol number of field columns, if genLayout is TRUE
#' @param serpentine a logical variable indicating whether plot number will be arranged as serpentine order
#' @param topToBottom a logical variable indicating whether plot number will be written from to to bottom
#' @param checkTrmtList NULL or a character vector indicating the names of the replicated treatment
#' @param testTrmtList NULL or a character vector indicating the names of the unreplicated treatment     
#' 
#' @return A list containing dataframe and statistical design array.
#'
#' @examples
#' randomizeAugmentedAL <- function(numCheck = 6, numTest = 12, numBlk = 3, numRep = 2,  
#'                                  genLayout = FALSE,
#'                                  checkTrmtList = paste("Check",1:6,sep=""), 
#'                                  testTrmtList = paste("Test",1:12,sep=""))
#'
# -------------------------------------------------------------------------------------

randomizeAugmentedAL <- function(numCheck, numTest, numBlk = 2, numRep = 2,  
                                 genLayout = FALSE, numRowPerBlk = 1, numColPerBlk,
                                 numRowPerRep = 1, numColPerRep, 
                                 numFieldRow = 1, numFieldCol, 
                                 serpentine = FALSE, topToBottom = TRUE, 
                                 checkTrmtList = NULL, testTrmtList = NULL) { #UseMethod("randomizeAugmentedAlpha") }
  
  # randomizeAugmentedAL.default <- function(numCheck, numTest, numBlk = 2, numRep = 2,  
  #                                          genLayout = FALSE, numRowPerBlk = 1, numColPerBlk,
  #                                          numRowPerRep = 1, numColPerRep, 
  #                                          numFieldRow = 1, numFieldCol, 
  #                                          serpentine = FALSE, topToBottom = TRUE, 
  #                                          checkTrmtList = NULL, testTrmtList = NULL) {
  
  numTrmt <- numTest + numCheck
  numPlots <- (numCheck * numRep) + numTest
  numPlotPerRep <- numCheck + (numTest/numRep) 
  blksize <- numPlotPerRep/numBlk
  numRepRow <- numFieldRow/numRowPerRep 
  numRepCol <- numFieldCol/numColPerRep                    
  numBlkRow <- numRowPerRep/numRowPerBlk                   
  numBlkCol <- numColPerRep/numColPerBlk
  
  tmpResult <- randomizeIBD(numTrmt = numCheck, 
                            numFieldRow = 1, numFieldCol = numCheck*numRep, 
                            numRep = numRep, numRowPerRep = 1, numColPerRep = numCheck,
                            numBlk = numBlk, blksize = numCheck/numBlk,
                            numRowPerBlk = 1, numColPerBlk = numCheck/numBlk,
                            trmtList = checkTrmtList)
  
  fieldbook1 <- tmpResult$fieldbook
  
  fieldbook1$Block <- rep(rep(c(1:numBlk), each = numCheck/numBlk),numRep)
  fieldbook1 <- fieldbook1[,c("ID", "ENTRY", "Block","REP")]
  fieldbook1$ENTRY <- as.numeric(fieldbook1$ENTRY) + numTest 
  
  fieldbook2 <- data.frame(PLOT = 1:numTest,
                           ENTRY = sample(1:numTest), 
                           Block = rep(1:numBlk, each = numTest/(numBlk*numRep)))
  
  fieldbook2 <- merge(data.frame(ID = testTrmtList, ENTRY = 1:numTest),
                      fieldbook2, by = "ENTRY")
  fieldbook2 <- fieldbook2[order(fieldbook2$PLOT),]
  row.names(fieldbook2) <- 1:numTest
  fieldbook2$REP <- rep(1:numRep, each = numTest/numRep)
  fieldbook2 <- fieldbook2[,c("ID", "ENTRY", "Block","REP")]
  
  fieldbook12 <- rbind(fieldbook1, fieldbook2)
  fieldbook12 <- fieldbook12[order(fieldbook12$REP, fieldbook12$Block),]
  
  fieldbook12$RAND <- as.vector(replicate(numBlk*numRep,sample(1:blksize)))
  fieldbook12 <- fieldbook12[order(fieldbook12$REP, fieldbook12$Block, fieldbook12$RAND),]
  row.names(fieldbook12) <- 1:numPlots
  fieldbook12$PlotNumPerRep <- rep(1:numPlotPerRep,numRep)
  fieldbook12$PlotNumber <- as.numeric(paste(fieldbook12$REP, paste(rep(0,nchar(nlevels(as.factor(fieldbook12$ID)))), collapse = ""), sep = "")) + fieldbook12$PlotNumPerRep
  fieldbook12 <- fieldbook12[,c("ENTRY", "PlotNumber", "Block", "REP", "ID")]
  
  if(!genLayout){
    tempLayout <- data.frame(index = rep(1:numPlotPerRep,numRep),fieldbook12)
    tempLayoutTrmt <- reshape(tempLayout[,c(1,5:6)], v.names = "ID", idvar = "index",
                              timevar = "REP", direction = "wide")
    tempLayoutTrmt <- tempLayoutTrmt[2:ncol(tempLayoutTrmt)]
    names(tempLayoutTrmt) <- paste("Rep",1:numRep, sep = "")
    
    tempLayoutBlock <- reshape(tempLayout[,c(1,4:5)], v.names = "Block", idvar = "index",
                               timevar = "REP", direction = "wide")
    tempLayoutBlock <- tempLayoutBlock[2:ncol(tempLayoutBlock)]
    names(tempLayoutBlock) <- paste("Rep",1:numRep, sep = "")
    
    trmtLayout <- list(tempLayoutTrmt, tempLayoutBlock)
    
    fieldbook <- fieldbook12[,c("ID", "ENTRY", "REP", "PlotNumber", "Block")]
    
  } else {
    fieldbook <- fieldbook12
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
        
        range <- getRange (index = k, topToBottom = topToBottom, serpentine = serpentine,
                           numRowPerGrp = numRowPerBlk, numColPerGrp = numColPerBlk,
                           numGrpRow = numBlkRow, numGrpCol = numBlkCol,
                           numRow = numRowPerRep, numCol = numColPerRep,
                           rowIndexLL = rowIndexLL, rowIndexUL = rowIndexUL,
                           colIndexLL = colIndexLL, colIndexUL = colIndexUL)
        
        rowIndexLL <- range$rowIndexLL
        rowIndexUL <- range$rowIndexUL
        colIndexLL <- range$colIndexLL
        colIndexUL <- range$colIndexUL
        
        if(k == 1){
          m <- ((j-1)*numPlotPerRep) + 1
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
    
    tempfbook <- as.data.frame.table(tempLayoutPlotNum)
    tempfbook$Var1 <- as.numeric(tempfbook$Var1)
    tempfbook$Var2 <- as.numeric(tempfbook$Var2)
    
    fieldbook <- merge(fieldbook, tempfbook, by.x = "PlotNumber", by.y = "Freq", all.x = TRUE)
    names(fieldbook)[6:7] <- c("ROW", "RANGE")
    
    trmtLayout <- list(tempLayoutTrmt, tempLayoutPlotNum, tempLayoutRep, tempLayoutBlock)
    
    fieldbook <- fieldbook[,c("ID", "ENTRY", "ROW", "RANGE", "REP", "PlotNumber", "Block")]
  }
  
  return(invisible(list(fieldbook = fieldbook, plan = trmtLayout)))
  
}

