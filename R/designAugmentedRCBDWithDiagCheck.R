# -------------------------------------------------------------------------------------
# Name             : designAugmentedRCBDWithDiagCheck
# Description      : Generate randomization for Augmented RCBD with Diagonal Checks
# R Version        : 4.1.2
# -------------------------------------------------------------------------------------
# Author           : Justine B. Dayrit
# Author Email     : j.bonifacio@irri.org
# Script Version   : 1.0.1
# Modified by      : Justine B. Dayrit
# Date Modified    : 2022.10.07
# Changes made     : updated matrix data to include trial i
# -------------------------------------------------------------------------------------
#' @title Augmented Randomized Complete Block Design with Diagonal Checks
#' @aliases designAugmentedRCBWithDiagCheck
#  @aliases designAugmentedRCBWithDiagCheck.default
#'  
#' @description Generate randomization and layout for Augmented Randomized Complete
#'       Block Design.
#'
#' @param numCheck number of replicated treatments
#' @param numTest number of unreplicated treatments
#' @param numSpatCheck number of spatial checks
#' @param trmtName NULL or a character string which indicate the name of the treatment to be displayed 
#'      in the output dataframe
#' @param numBlk number of blocks
#' @param numTrial number of trials (randomization set-ups with the complete set of entries)
#' @param numRowPerBlk number of rows per block, if genLayout is TRUE 
#' @param numFieldRow number of field rows, if genLayout is TRUE
#' @param serpentine a logical variable indicating whether plot number will be arranged as serpentine order
#' @param topToBottom a logical variable indicating whether plot number will be written from to to bottom
#' @param checkTrmtList NULL or a character vector indicating the names of the replicated treatment
#' @param testTrmtList NULL or a character vector indicating the names of the unreplicated treatment     
#' @param spatCheckTrmtList NULL or a character vector indicating the names of the spatial checks     
#' @param display a logical variable indicating whether randomization parameters will be displayed
#' 
#' @details 
#' The parameters \code{numFieldRow} and \code{numRowPerBlk} must be specified.  
#' Values of \code{numRowPerBlk} should be a factor of \code{numFieldRow}.
#' 
#' If \code{checkTrmtList} is a character vector, the length should be equal to the \code{numCheck}.
#' If \code{testTrmtList} is a character vector, the length should be equal to the \code{numTest}.
#' If \code{spatCheckTrmtList} is a character vector, the length should be equal to the \code{numSpatCheck}.
#' 
#' @return A list containing the following components:
#' \item{fieldbook}{a data frame} 
#' \item{plan}{a list containing the following components:}
#' \item{TrmtLayout}{a list whose length is equal to the number of trials containing the treatment layout for each trial}
#' \item{PlotNumLayout}{a matrix containing the plot number information of the experiment}
#' \item{BlockLayout}{a matrix containing the replication information of the experiment}
#' 
#' @examples
#' ## Generate randomization of 6 entries replicated 4 times, 8 unreplicated entries and 2 spatial checks
#' ## using Augmented Randomized Complete Block Design with Diagonal checks
# augRCBDd1 <- designAugmentedRCBWithDiagCheck(numCheck = 6, numTest = 8, numSpatCheck =2,
#                                              trmtName = "Variety", numBlk = 4,
#                                              numRowPerBlk = 1, numFieldRow = 4,
#                                              serpentine = TRUE, topToBottom = FALSE,
#                                              numTrial = 2)
#'                                
# -------------------------------------------------------------------------------------

#install.packages("dplyr")
#library(dplyr)
#library(crayon)

designAugmentedRCBWithDiagCheck <- function(numCheck, numTest, numSpatCheck, trmtName = NULL, numBlk = 2, numTrial = 1, 
                                            numRowPerBlk = 1, numFieldRow = 1, serpentine = FALSE, topToBottom = TRUE, 
                                            checkTrmtList = NULL, testTrmtList = NULL, 
                                            spatCheckTrmtList = NULL, display = TRUE) { #UseMethod("designAugmentedRCB") }

#designAugmentedRCBWithDiagCheck.default <- function(numCheck, numTest, numSpatCheck, trmtName = NULL, numBlk = 2, numTrial = 1, 
#                                                    numRowPerBlk = 1, numFieldRow = 1, 
#                                                    serpentine = FALSE, topToBottom = TRUE, checkTrmtList = NULL, testTrmtList = NULL, 
#                                                    spatCheckTrmtList = NULL, display = TRUE) {
  
  # check user input
  if (numBlk < 2) { 
    stop("ERROR: The number of blocks should be greater than or equal to 2.") }  # -- check the number of replicates
  
  if (numSpatCheck < 1) { 
    stop("ERROR: There should be at least 1 spatial check.")}
  
  if (numSpatCheck > numBlk) { 
    stop("ERROR: The number of spatial checks should be less than the total number of blocks.") }
  
  # -- determine if the number of elements in checkTrmtList/testTrmtList is equal to numCheck/numTest
  if (is.null(checkTrmtList)) {
    checkTrmtList <- paste("Check", 1:numCheck, sep = "")
    checkList <- FALSE
  } else {
    if (length(checkTrmtList) != numCheck) { 
      stop("ERROR: The number of elements of the arg 'checkTrmtList' is not equal to numCheck.") }
    checkList <- TRUE
  }
  
  if (is.null(testTrmtList)) { 
    testTrmtList <- paste("Test",1:numTest, sep = "")
    testList <- FALSE
  } else {
    if (length(testTrmtList) != numTest) { 
      stop("ERROR: The number of elements of the arg 'testTrmtList' is not equal to numTest.") }
    testList <- TRUE
  }
  
  if (is.null(spatCheckTrmtList)) {
    spatCheckTrmtList <- paste("SpatialCheck", 1:numSpatCheck, sep = "")
    spatCheckList <- FALSE
  } else {
    if (length(spatCheckTrmtList) != numSpatCheck) { 
      stop("ERROR: The number of elements of the arg 'spatCheckTrmtList' is not equal to numSpatCheck") }
    spatCheckList <- TRUE
  }

  if (numTest%%numBlk != 0) { 
    stop("ERROR: The number of test treatments should be divisible by the number of blocks.")  }
  
  numTestPerBlk <- numTest/numBlk                    # -- determine the number of test entries per blk
  numTrmt <- numCheck + numTest                      # -- determine the total number of treatment in the experiment
  numPlots <- (numCheck * numBlk) + numTest          # -- determine the total number of experimental units 
  errorDF <- (numBlk - 1) * (length(checkTrmtList) - 1)  # -- compute the error df for a complete experiment
  numPlotPerBlk <- numPlots/numBlk                   # -- assume equal number of plots per block
  
  numFieldCol <- numPlots/numFieldRow                  # -- determine the number of columns in the design
  numColPerBlk <- numPlotPerBlk/numRowPerBlk           # -- determine the number of columns per block
  
  # -- determine the orientation of blocks
  numBlkRow <- numFieldRow/numRowPerBlk
  numBlkCol <- numFieldCol/numColPerBlk
  
  if (numBlkRow * numBlkCol != numBlk) { 
    stop("ERROR: The number of field rows should be divisible by the number of rows per block.") }
  
  #if (numRowPerBlk == 1 || numRowPerBlk == numPlotPerBlk) serpentine <- FALSE
  
  if (numFieldRow%%numRowPerBlk != 0) { 
    stop("ERROR: The number of field rows should be divisible by the number of rows per block.") }
  
  # generate randomization and fieldbook
  # capture.output(result <- designRCBD(generate = list(EntryNo = 1:numPlotPerBlk), r, trial, numFieldRow, rowPerBlk, serpentine, display = FALSE))
  
  capture.output(tmpResult <- designRCBD(generate = list(Entry = 1:numPlotPerBlk), 
                                         numBlk = numBlk, numTrial = numTrial,
                                         genLayout = TRUE, 
                                         numFieldRow = numFieldRow, 
                                         numRowPerBlk = numRowPerBlk,
                                         serpentine = FALSE,
                                         topToBottom = topToBottom))
  
  fbook <- tmpResult$fieldbook[order(tmpResult$fieldbook$Trial, 
                                     tmpResult$fieldbook$PlotNumber),]
  
  #-- recoding the entry number for augmented design in RCB
  fbook[,"Entry"] <- as.numeric(fbook[,"Entry"])
  
  newfbook <- NULL
  newfbook1 <- NULL
  
  tmpCheckIDNum <- sample(numPlotPerBlk, numCheck)
  tmpTestIDNum <- setdiff(1:numPlotPerBlk, tmpCheckIDNum)
  
  for (i in 1:numTrial) {
    tmpTestList <- sample(testTrmtList, numTest)
    
    from <- c(tmpCheckIDNum, tmpTestIDNum)
    to <- c(checkTrmtList, rep(NA, each = numTestPerBlk))
    
    tmp1fbook <- fbook[fbook$Trial == i,]
    tmp1fbook$tmpEntryCode <- recode(tmp1fbook$Entry, !!! setNames(to, from))
    tmp1fbook[is.na(tmp1fbook$tmpEntryCode),"tmpEntryCode"] <- tmpTestList
    tmp1fbook$Entry <- tmp1fbook$tmpEntryCode
    tmp1fbook$tmpEntryCode <- NULL
    
    if(i == 1){
      if(numRowPerBlk == 1){
        numBlkPerRow <- numFieldCol/numColPerBlk
        numFieldCol <- numFieldCol + numBlkPerRow
        numColPerBlk <- numColPerBlk + 1
        percentCheck <- (numBlkPerRow*numFieldRow)/(numFieldCol*numFieldRow)
      } else if(numRowPerBlk == numPlotPerBlk){
        numBlkPerCol <- numFieldRow/numRowPerBlk
        numFieldRow <- numFieldRow + numBlkPerCol
        numRowPerBlk <- numRowPerBlk + 1
        percentCheck <- (numBlkPerCol*numFieldCol)/(numFieldRow*numFieldCol)
      } else {
        stop(paste("This feature is not yet available"))
      }
      
      numBlkRow <- numFieldRow/numRowPerBlk 
      numBlkCol <- numFieldCol/numColPerBlk
      numPlotPerBlk <- numPlotPerBlk + 1
    }
    
    if(numRowPerBlk == 1){
      dlayout <- diag.layout(numFieldCol, numFieldRow, percentCheck, plot = FALSE)
      temptrmt <- as.vector(t(dlayout))
    } else if(numRowPerBlk == numPlotPerBlk){
      dlayout <- diag.layout(numFieldRow, numFieldCol, percentCheck, plot = FALSE)
      dlayout <- t(dlayout)
      temptrmt <- as.vector(dlayout)
    }
    
    if((sum(dlayout)%%numSpatCheck) == 0){
      tempSpat <- rep(spatCheckTrmtList,each = (sum(dlayout)/numSpatCheck))
      tempSpatOrder <- data.frame(tempSpat , tempPlotNum = sample(length(tempSpat), length(tempSpat), replace = FALSE))
      tempSpatOrder <- tempSpatOrder[order(tempSpatOrder[,"tempPlotNum"]),]
    } else{
      tempSpat <- c(rep(spatCheckTrmtList,each = floor(sum(dlayout)/numSpatCheck)),sample(spatCheckTrmtList, (sum(dlayout)%%numSpatCheck), replace = FALSE))
      tempSpatOrder <- data.frame(tempSpat , tempPlotNum = sample(length(tempSpat), length(tempSpat), replace = FALSE))
      tempSpatOrder <- tempSpatOrder[order(tempSpatOrder[,"tempPlotNum"]),]
    }
    
    temptrmt[temptrmt == 1] <- tempSpatOrder$tempSpat
    #temptrmt[temptrmt == 1] <- as.numeric(as.vector(sample(as.factor(spatCheckNumS:spatCheckNumE), sum(dlayout), replace = TRUE)))
    temptrmt[temptrmt == 0] <- as.vector(tmp1fbook[1:nrow(tmp1fbook),"Entry"])
    
    if(numRowPerBlk == 1){
      tempnewfbook1 <- data.frame("Trial" = i, "Block" = rep(1:numBlk,each = numPlotPerBlk), "Entry" = temptrmt, 
                              "PlotNumber" = 1:length(dlayout), "FieldRow" = rep(1:numFieldRow, each = numFieldCol),
                              "FieldColumn" = rep(1:numFieldCol, numFieldRow))
    } else if(numRowPerBlk == numPlotPerBlk){
      tempnewfbook1 <- data.frame("Trial" = i, "Block" = rep(1:numBlk,each = numPlotPerBlk), "Entry" = temptrmt, 
                              "PlotNumber" = 1:length(dlayout), "FieldRow" = rep(1:numFieldRow, numFieldCol),
                              "FieldColumn" = rep(1:numFieldCol, each = numFieldRow))
    }
    
    tempnewfbook1$PlotNumber <- as.numeric(paste(tempnewfbook1$Block, paste(rep(0,nchar(length(dlayout)/numBlk)), collapse=""), sep = "")) + 1:(length(dlayout)/numBlk)
    
    if(numRowPerBlk == 1){
      if(!topToBottom & serpentine){
        tempPlot <- matrix(tempnewfbook1$PlotNumber, nrow = numFieldRow, ncol = numFieldCol, byrow = TRUE)
        tempBlk <- matrix(tempnewfbook1$Block, nrow = numFieldRow, ncol = numFieldCol, byrow = TRUE)
        for (z in seq(2, numFieldRow, by = 2)) { 
          tempPlot[z,] <- rev(tempPlot[z,]) 
          tempBlk[z,] <- rev(tempBlk[z,]) }
        tempnewfbook1$PlotNumber <- as.vector(t(tempPlot))
        tempnewfbook1$Block <- as.vector(t(tempBlk))
      }
    } else if(numRowPerBlk == numPlotPerBlk){
      if (topToBottom & serpentine){
        tempPlot <- matrix(tempnewfbook1$PlotNumber, nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)
        tempBlk <- matrix(tempnewfbook1$Block, nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)
        for (z in seq(2, numFieldCol, by = 2)) { 
          tempPlot[,z] <- rev(tempPlot[,z]) 
          tempBlk[,z] <- rev(tempBlk[,z]) }
        tempnewfbook1$PlotNumber <- as.vector(tempPlot)
        tempnewfbook1$Block <- as.vector(tempBlk)
      }
    }
    
    newfbook1 <- rbind(newfbook1, tempnewfbook1)
    
    if(topToBottom){
      if(numRowPerBlk == numPlotPerBlk){
        tempLayoutTrmt <- matrix(newfbook1[newfbook1$Trial == i, "Entry"], nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)  
      } else{
        tempLayoutTrmt <- matrix(newfbook1[newfbook1$Trial == i, "Entry"], nrow = numFieldRow, ncol = numFieldCol, byrow = TRUE)  
      }
    } else{
      if(numRowPerBlk == 1){
        tempLayoutTrmt <- matrix(newfbook1[newfbook1$Trial == i, "Entry"], nrow = numFieldRow, ncol = numFieldCol, byrow = TRUE) 
      } else{
        tempLayoutTrmt <- matrix(newfbook1[newfbook1$Trial == i, "Entry"], nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)  
      }
    }
    
    tmpResult$plan$TrmtLayout[[i]] <- tempLayoutTrmt
    dimnames(tmpResult$plan$TrmtLayout[[i]]) <- list(paste("FieldRow", 1:numFieldRow,sep = ""),
                                                     paste("FieldColumn", 1:numFieldCol, sep = ""))
    if(i == 1){
      if(topToBottom){
        if(numRowPerBlk == numPlotPerBlk){
          tmpResult$plan$PlotNumLayout <- matrix(newfbook1$PlotNumber, nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)
          tmpResult$plan$BlockLayout <-  matrix(newfbook1$Block, nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)  
        } else{
          tmpResult$plan$PlotNumLayout <- matrix(newfbook1$PlotNumber, nrow = numFieldRow, ncol = numFieldCol, byrow = TRUE)
          tmpResult$plan$BlockLayout <-  matrix(newfbook1$Block, nrow = numFieldRow, ncol = numFieldCol, byrow = TRUE)
        }
      } else{
        if(numRowPerBlk == 1){
          tmpResult$plan$PlotNumLayout <- matrix(newfbook1$PlotNumber, nrow = numFieldRow, ncol = numFieldCol, byrow = TRUE)
          tmpResult$plan$BlockLayout <-  matrix(newfbook1$Block, nrow = numFieldRow, ncol = numFieldCol, byrow = TRUE)
        } else{
          tmpResult$plan$PlotNumLayout <- matrix(newfbook1$PlotNumber, nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)
          tmpResult$plan$BlockLayout <-  matrix(newfbook1$Block, nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)
        }
      }
      
      dimnames(tmpResult$plan$PlotNumLayout) <- list(paste("FieldRow", 1:numFieldRow,sep = ""),
                                                     paste("FieldColumn", 1:numFieldCol, sep = ""))
      dimnames(tmpResult$plan$BlockLayout) <- list(paste("FieldRow", 1:numFieldRow,sep = ""),
                                                   paste("FieldColumn", 1:numFieldCol, sep = ""))
    }
  }
  
  newfbook1 <- newfbook1[order(newfbook1$Trial, newfbook1$PlotNumber),]
  
  newfbook1 <- as.data.frame(newfbook1 %>% 
                         group_by_at(vars(Entry,Trial)) %>% 
                         mutate(Rep=1:n(),.after = 1))
  
  if (display) {
    cat(toupper("Design Properties:"), "\n", sep = "")
    cat("\t", "Augmented Randomized Complete Block Design with Diagonal Checks", "\n\n", sep = "")
    cat(toupper("Design Parameters"), "\n", sep = "")
    cat("\t", "Number of Trials = ", numTrial, "\n", sep = "")
    cat("\t", "Number of Replicated Treatments = ", numCheck, "\n", sep = "")
    cat("\t", "Levels of Replicated Treatments = ", sep = "")
    
    if (numCheck <= 5) { 
      cat(paste(checkTrmtList, collapse = ", ", sep = ""), "\n", sep = "") 
    } else {
      cat(paste(checkTrmtList[1:3], collapse = ", ", sep = ""), sep = "") 
      cat(", ..., ", checkTrmtList[numCheck], "\n", sep = "") 
    }
    
    cat("\t", "Number of Blocks = ", numBlk, "\n", sep = "")
    cat("\t", "Number of UnReplicated Treatments = ", numTest, "\n", sep = "")
    cat("\t", "Levels of UnReplicated Treatments = ", sep = "")
    
    if (numTest <= 5) { 
      cat(paste(testTrmtList, collapse = ", ", sep = ""), "\n", sep = "") 
    } else {
      cat(paste(testTrmtList[1:3], collapse = ", ", sep = ""), sep = "") 
      cat(", ..., ", testTrmtList[numTest], "\n", sep = "") 
    }
    
    cat("\t", "Number of Spatial Checks = ", numSpatCheck, "\n\n", sep = "")
    
    cat("\t", "Number of Field Rows = ", numFieldRow,"\n", sep = "")
    cat("\t", "Number of Field Columns = ", numFieldCol, "\n\n", sep = "")
    
  } ## end stmt -- if (display)
  
  if (errorDF < 12) { cat(green("WARNING: Too few error df.","\n\n")) }
  if (numPlotPerBlk > 30) { cat(green("WARNING: The block size is too large. You might want to consider increasing the number of blocks.","\n\n"))}
    
  if (testList == FALSE & checkList == FALSE){
    if (!is.null(trmtName)) {
      names(newfbook1)[match("Entry", names(newfbook1))] <- trmtName
    }
  } else {
    entryList <- data.frame(Entry = c(testTrmtList,checkTrmtList,spatCheckTrmtList), EntryType = c(rep("Test",numTest),rep("Check",numCheck),rep("Spatial Check",numSpatCheck)))
    newfbook1 <- merge(newfbook1,entryList)
    newfbook1 <- newfbook1[,c("Trial","Rep","Block","Entry","EntryType","PlotNumber","FieldRow","FieldColumn")]
    newfbook1 <- newfbook1[order(newfbook1$Trial, newfbook1$PlotNumber), ]
    
    if (!is.null(trmtName)) {
      names(newfbook1)[match("Entry", names(newfbook1))] <- trmtName
      names(newfbook1)[match("EntryType", names(newfbook1))] <- paste(trmtName,"Type",sep="")}
  }
  
  return(invisible(list(fieldbook = newfbook1, plan = tmpResult$plan)))
  
} ## end stmt -- designAugmentedRCBWithDiagCheck function

