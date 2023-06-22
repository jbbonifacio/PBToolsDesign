# -------------------------------------------------------------------------------------
# Name             : designAugmentedRCBD
# Description      : Generate randomization for Augmented RCBD
# R Version        : 4.0.1
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 4.1.1
# Modified by      : Justine B. Bonifacio
# Date Modified    : 2020.03.26
# Changes made     : add entry type when checkTrmtList and testTrmtList are not NULL
#                    instead of using "new" for identifying test treatments, use "Test"
# -------------------------------------------------------------------------------------
#' @title Augmented Randomized Complete Block Design
#' @aliases designAugmentedRCB
#  @aliases designAugmentedRCB.default
#'  
#' @description Generate randomization and layout for Augmented Randomized Complete
#'       Block Design.
#'
#' @param numCheck number of replicated treatments
#' @param numTest number of unreplicated treatments
#' @param trmtName NULL or a character string which indicate the name of the treatment to be displayed 
#'      in the output dataframe
#' @param numBlk number of blocks
#' @param numTrial number of trials (randomization set-ups with the complete set of entries)
#' @param genLayout logical, whether a layout of the design will be generated or not
#' @param numRowPerBlk number of rows per block, if genLayout is TRUE 
#' @param numFieldRow number of field rows, if genLayout is TRUE
#' @param serpentine a logical variable indicating whether plot number will be arranged as serpentine order
#' @param topToBottom a logical variable indicating whether plot number will be written from to to bottom
#' @param checkTrmtList NULL or a character vector indicating the names of the replicated treatment
#' @param testTrmtList NULL or a character vector indicating the names of the unreplicated treatment     
#' @param display a logical variable indicating whether randomization parameters will be displayed
#' 
#' @details 
#' If \code{genLayout} is \code{TRUE}, then parameters \code{numFieldRow} and \code{numRowPerBlk} must be specified.  
#' Values of \code{numRowPerBlk} should be a factor of \code{numFieldRow}.
#' 
#' If \code{checkTrmtList} is a character vector, the length should be equal to the \code{numCheck}.
#' If \code{testTrmtList} is a character vector, the length should be equal to the \code{numTest}.
#' 
#' @return A list containing the following components:
#' \item{fieldbook}{a data frame} 
#' \item{plan}{a data frame, if \code{genLayout} is \code{FALSE} or a list containing the following components 
#' if \code{genLayout} is \code{TRUE}:}
#' \item{TrmtLayout}{a list whose length is equal to the number of trials containing the treatment layout for each trial}
#' \item{PlotNumLayout}{a matrix containing the plot number information of the experiment}
#' \item{BlockLayout}{a matrix containing the replication information of the experiment}
#' 
#' @examples
#' ## Generate randomization of 6 entries replicated 4 times and 8 unreplicated entries using 
#' ## Augmented Randomized Complete Block Design
#' augRCBD1 <- designAugmentedRCB(numCheck = 6, numTest = 8, trmtName = "Variety", numBlk = 4,
#'     numTrial = 2, genLayout = FALSE)
#'            
#' augRCBD2 <- designAugmentedRCB(numCheck = 6, numTest = 8, trmtName = "Variety", numBlk = 4,
#'     numTrial = 2, genLayout = TRUE, numRowPerBlk = 2, numFieldRow = 4, serpentine = TRUE)
#'                                
# -------------------------------------------------------------------------------------

#install.packages("dplyr")
#library(dplyr)
#library(crayon)

designAugmentedRCB <- function(numCheck, numTest, trmtName = NULL, numBlk = 2, numTrial = 1, 
                               genLayout = FALSE, numRowPerBlk = 1, numFieldRow = 1, 
                               serpentine = FALSE, topToBottom = TRUE, checkTrmtList = NULL, testTrmtList = NULL, 
                               display = TRUE) { #UseMethod("designAugmentedRCB") }

#designAugmentedRCB.default <- function(numCheck, numTest, trmtName = NULL, numBlk = 2, numTrial = 1, 
#                                       genLayout = FALSE, numRowPerBlk = 1, numFieldRow = 1, 
#                                       serpentine = FALSE, topToBottom = TRUE, checkTrmtList = NULL, testTrmtList = NULL, 
#                                       display = TRUE) {
  
  # check user input
  if (numBlk < 2) { 
    stop("ERROR: The number of blocks should be greater than or equal to 2.") }  # -- check the number of replicates
  
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

  if (numTest%%numBlk != 0) { 
    stop("ERROR: The number of test treatments should be divisible by the number of blocks.")  }
  
  numTestPerBlk <- numTest/numBlk                    # -- determine the number of test entries per blk
  numTrmt <- numCheck + numTest                      # -- determine the total number of treatment in the experiment
  numPlots <- (numCheck * numBlk) + numTest          # -- determine the total number of experimental units 
  errorDF <- (numBlk - 1) * (length(checkTrmtList) - 1)  # -- compute the error df for a complete experiment
  numPlotPerBlk <- numPlots/numBlk                   # -- assume equal number of plots per block
  
  if (genLayout) {
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
  }
  
  # generate randomization and fieldbook
  # capture.output(result <- designRCBD(generate = list(EntryNo = 1:numPlotPerBlk), r, trial, numFieldRow, rowPerBlk, serpentine, display = FALSE))
  
  capture.output(tmpResult <- designRCBD(generate = list(Entry = 1:numPlotPerBlk), 
                                         numBlk = numBlk, numTrial = numTrial,
                                         genLayout = genLayout, 
                                         numFieldRow = numFieldRow, 
                                         numRowPerBlk = numRowPerBlk,
                                         serpentine = serpentine,
                                         topToBottom = topToBottom))
  
  if(genLayout) {
    fbook <- tmpResult$fieldbook[order(tmpResult$fieldbook$Trial,
                                       tmpResult$fieldbook$FieldRow,
                                       tmpResult$fieldbook$FieldColumn),]
  } else {
    fbook <- tmpResult$fieldbook[order(tmpResult$fieldbook$Trial,
                                       tmpResult$fieldbook$PlotNumber),]
  }

  #-- recoding the entry number for augmented design in RCB
  fbook[,"Entry"] <- as.numeric(fbook[,"Entry"])
  
  newfbook <- NULL
  
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
    newfbook <- rbind(newfbook,tmp1fbook)
    
    if(genLayout){
      tmp1 <- as.vector(tmpResult$plan$TrmtLayout[[i]])
      tmp2 <- recode(tmp1, !!! setNames(to, from))
      tmp2[is.na(tmp2)] <- tmpTestList
      tmpResult$plan$TrmtLayout[[i]] <- matrix(tmp2, nrow = numFieldRow, ncol = numFieldCol, byrow = FALSE)
      dimnames(tmpResult$plan$TrmtLayout[[i]]) <- list(paste("FieldRow", 1:numFieldRow,sep = ""),
                                                       paste("FieldColumn", 1:numFieldCol, sep = ""))
    }
  }
  
  newfbook <- newfbook[order(fbook$Trial, fbook$PlotNumber),]
  
  fbook <- as.data.frame(newfbook %>% 
                         group_by_at(vars(Entry,Trial)) %>% 
                         mutate(Rep=1:n(),.after = 1))
  
  if (display) {
    cat(toupper("Design Properties:"), "\n", sep = "")
    cat("\t", "Augmented Randomized Complete Block Design (Augmented RCBD)", "\n\n", sep = "")
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
      cat(", ..., ", testTrmtList[numTest], "\n\n", sep = "") 
    }
    
    if (genLayout) {
      cat("\t", "Number of Field Rows = ", numFieldRow,"\n", sep = "")
      cat("\t", "Number of Field Columns = ", numFieldCol, "\n\n", sep = "")
    }
    
  } ## end stmt -- if (display)
  
  if (errorDF < 12) { cat(green("WARNING: Too few error df.","\n\n")) }
  if (numPlotPerBlk > 30) { cat(green("WARNING: The block size is too large. You might want to consider increasing the number of blocks.","\n\n"))}
  
  if (!genLayout) {
    tmp1 <- data.frame(index = rep(1:numPlotPerBlk, numBlk*numTrial), fbook)
    plan <- reshape(tmp1, v.names = "Entry", idvar = c("index","Trial"),
                    timevar = "Block", direction = "wide",
                    drop = c("PlotNumber", "Rep"))
    plan <- plan[2:ncol(plan)]
    names(plan)[2:ncol(plan)] <- paste("Block",1:numBlk, sep = "")
    rownames(plan) <- 1:nrow(plan)
    
    if (testList == FALSE & checkList == FALSE){
      if (!is.null(trmtName)) {
        names(fbook)[match("Entry", names(fbook))] <- trmtName
      }
    } else {
      entryList <- data.frame(Entry = c(testTrmtList,checkTrmtList), EntryType = c(rep("Test",numTest),rep("Check",numCheck)))
      fbook <- merge(fbook,entryList)
      fbook <- fbook[,c("Trial","Rep","Block","Entry","EntryType","PlotNumber")]
      fbook <- fbook[order(fbook$Trial, fbook$PlotNumber), ]
      
      if (!is.null(trmtName)) {
        names(fbook)[match("Entry", names(fbook))] <- trmtName
        names(fbook)[match("EntryType", names(fbook))] <- paste(trmtName,"Type",sep="")}
    }
    
    return(invisible(list(fieldbook = fbook, plan = plan)))
  } else {
    
    if (testList == FALSE & checkList == FALSE){
      if (!is.null(trmtName)) {
        names(fbook)[match("Entry", names(fbook))] <- trmtName
      }
    } else {
      entryList <- data.frame(Entry = c(testTrmtList,checkTrmtList), EntryType = c(rep("Test",numTest),rep("Check",numCheck)))
      fbook <- merge(fbook,entryList)
      fbook <- fbook[,c("Trial","Rep","Block","Entry","EntryType","PlotNumber","FieldRow","FieldColumn")]
      fbook <- fbook[order(fbook$Trial, fbook$PlotNumber), ]
      
      if (!is.null(trmtName)) {
        names(fbook)[match("Entry", names(fbook))] <- trmtName
        names(fbook)[match("EntryType", names(fbook))] <- paste(trmtName,"Type",sep="")}
    }
    
    return(invisible(list(fieldbook = fbook, plan = tmpResult$plan)))
  }
  
} ## end stmt -- designAugmentedRCBD function

