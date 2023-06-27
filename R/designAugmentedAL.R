# -------------------------------------------------------------------------------------
# Name             : designAugmentedAL
# Description      : Generate randomization for Augmented Alpha Lattice
# R Version        : 4.1.0
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles
# Author Email     : a.gulles@irri.org
# Script Version   : 1.1.0
# Modified by      : Justine B. Bonifacio
# Date Modified    : 2021.09.09
# Changes made     : changed the implementation of randomization, 3 steps
#                    step 1: randomize the checks, using ibDiGGer
#                    step 2: randomize the test, using sample
#                    step 3: randomize each block, using sample
# -------------------------------------------------------------------------------------
#' @title Augmented Alpha Lattice Design
#' @aliases designAugmentedAL
#  @aliases designAugmentedAL.default
#'  
#' @description Generate randomization and layout for Augmented Alpha Lattice Design.
#'
#' @param numCheck number of replicated treatments
#' @param numTest number of unreplicated treatments
#' @param trmtName NULL or a character string which indicate the name of the treatment to be displayed in the output dataframe
#' @param numBlk number of blocks per replicate
#' @param numRep number of replicates
#' @param numTrial number of trials (randomization set-ups with the complete set of entries)
#' @param genLayout logical, whether a layout of the design will be generated or not
#' @param numRowPerBlk number of rows per block, if genLayout is TRUE 
#' @param numRowPerRep number of rows per replicate
#' @param numFieldRow number of field rows, if genLayout is TRUE
#' @param serpentine a logical variable indicating whether plot number will be arranged as serpentine order
#' @param topToBottom a logical variable indicating whether plot number will be written from to to bottom
#' @param checkTrmtList NULL or a character vector indicating the names of the replicated treatment
#' @param testTrmtList NULL or a character vector indicating the names of the unreplicated treatment     
#' @param display a logical variable indicating whether randomization parameters will be displayed
#' 
#' @details 
#' If \code{genLayout} is \code{TRUE}, then parameters \code{numFieldRow}, \code{numRowPerRep}  and \code{numRowPerBlk} 
#' must be specified.  Values of \code{numRowPerBlk} should be a factor of \code{numRowPerRep} while values of 
#' \code{numRowPerRep} should be a factor of \code{numFieldRow}.
#' 
#' If \code{checkTrmtList} is a character vector, the length should be equal to the \code{numCheck}.
#' If \code{testTrmtList} is a character vector, the length should be equal to the \code{numTest}.
#' 
#' @return A list containing the following components:
#' \item{fieldbook}{a data frame}
#' \item{plan}{a list containing the following components: TrmtLayout and BlockLayout
#' if \code{genLayout} is \code{FALSE} or a list containing the following components
#' if \code{genLayout} is \code{TRUE}:}
#' \item{TrmtLayout}{a data frame containing the treatment layout of the experiment, 
#'      if \code{genLayout} is \code{FALSE} or a list whose length is equal to the 
#'      number of trials containing the treatment layout for each trial, 
#'      if \code{genLayout} is \code{TRUE}} 
#' \item{PlotNumLayout}{a matrix containing the plot number layout of the experiment} 
#' \item{RepLayout}{a matrix containing the replication layout of the experiment} 
#' \item{BlockLayout}{a data frame containing the block layout of the experiment,
#'      if \code{genLayout} is \code{FALSE} or a matrix containing the block layout of the experiment,
#'      if \code{genLayout} is \code{TRUE}}
#' 
#' @examples
#' ## Generate randomization of 8 entries replicated 2 times and 16 unreplicated entries using 
#' ## Augmented Alpha Lattice Design
#' augAL1 <- designAugmentedAL(numCheck = 8, numTest = 16, numBlk = 4, numRep = 2,  
#'                             numTrial = 1, genLayout = FALSE)
#'            
#' augAL2 <- designAugmentedAL(numCheck = 8, numTest = 16, trmtName = "Variety", numBlk = 4, 
#'                             numRep = 2, numTrial = 2, genLayout = TRUE, numRowPerBlk = 4, 
#'                             numRowPerRep = 4, numFieldRow = 8, 
#'                             serpentine = TRUE, topToBottom = TRUE)
#'                                                                
# -------------------------------------------------------------------------------------

designAugmentedAL <- function(numCheck, numTest, trmtName = NULL, numBlk = 2, numRep = 2,  
                              numTrial = 1, genLayout = FALSE, numRowPerBlk = 1, numRowPerRep = 1, numFieldRow = 1, 
                              serpentine = FALSE, topToBottom = TRUE, checkTrmtList = NULL, testTrmtList = NULL, 
                              display = TRUE) { #UseMethod("designAugmentedRCB") }
        
# designAugmentedAL.default <- function(numCheck, numTest, trmtName = NULL, numBlk = 2,  numRep = 2, 
#                                       numTrial = 1, genLayout = FALSE, numRowPerBlk = 1, numRowPerRep = 1, numFieldRow = 1, 
#                                       serpentine = FALSE, topToBottom = TRUE, checkTrmtList = NULL, testTrmtList = NULL, 
#                                       display = TRUE) { 

     # Check input of the user:
        if (numRep < 2) { stop("ERROR: The number of replicates should be greater than or equal to 2.")}
        
        if (numTest%%numRep != 0) { stop("ERROR: The number of test entries should be divisible by the number of replicates.")}
        
        if (numBlk == 1) { stop("ERROR: The number of blocks per replicate should be greater than 1.") }
        
        flag1 <- primeNumber(numTest)
        if(flag1 == 1) { stop("ERROR: The number of test entries should not be a prime number.")}
        
        flag2 <- primeNumber(numCheck)
        if(flag2 == 1) { stop("ERROR: The number of check entries should not be a prime number.")}
        
        if(numTest < 4){ stop("ERROR: The number of test entries should be greater than or equal to 4.")}
        if(numCheck < 4){ stop("ERROR: The number of check entries should be greater than or equal to 4.")}
        
        if (numCheck/numBlk < 2){ stop("ERROR: The number of check entries per block should be greater than or equal to 2.")}
        if ((numTest/numRep)%%numBlk != 0) { stop("ERROR: The number of test entries should be divisible by the number of replicates and number of blocks.")}
        
        treatList <- TRUE
        # --- determine if the number of elements in checkTrmt/newTrmt is equal to numCheck/numNew
        if (is.null(checkTrmtList)) { 
          checkTrmtList <- paste("Check",1:numCheck, sep = "")
          treatList <- FALSE
        } else {
                if (length(checkTrmtList) != numCheck) { stop("ERROR: The number of elements of the arg 'checkTrmtList' is not equal to numCheck.") }          
        }
        if (is.null(testTrmtList)) { 
          testTrmtList <- paste("Test",1:numTest, sep = "")
          treatList <- FALSE
        } else {
                if (length(testTrmtList) != numTest) { stop("ERROR: The number of elements of the arg 'testTrmtList' is not equal to numTest.") }          
        }
        
        checkTrmtList <- as.character(checkTrmtList)
        testTrmtList <- as.character(testTrmtList)
        
        # determine the total number of experimental units    
        numTrmt <- numCheck + numTest 
        numPlots <- (numCheck * numRep) + numTest
        
        if (!genLayout) {
                numFieldRow <- 1
                numRowPerRep <- 1
                numRowPerBlk <- 1
        }
        
        #if (numPlots > 1500) { stop("ERROR: The maximum number of experimental units that can be generated is 1500.") }
        if (numPlots %% numFieldRow != 0) { stop("ERROR: Total number of experimental units should be divisible by the number of field rows.") }
        numFieldCol <- numPlots/numFieldRow
        
        if (numPlots%%numRep != 0) { stop("ERROR: The total number of experimental units should be divisible by the number of replicates.") }
        numPlotPerRep <- numPlots/numRep
        
        if (numPlotPerRep%%numRowPerRep != 0) { stop("ERROR: The total number of experimental units per replicated should be divisible by the number of rows per replicates.") }
        numColPerRep <- numPlotPerRep/numRowPerRep
        
        #determine block size
        blksize <- numPlotPerRep/numBlk
        
        if (blksize%%numRowPerBlk != 0) { stop("ERROR: The total number of plots per block should be divisible by the number of rows per block.") }
        numColPerBlk <- blksize/numRowPerBlk
        
        fieldbook <- NULL
        trmtLayout <- NULL
        blkLayout <- NULL
        plotNumLayout <- NULL
        repLayout <- NULL
        tmpResult2 <- NULL
        
        for (i in 1:numTrial) {
                tmpResult <- randomizeAugmentedAL(numCheck = numCheck, numTest = numTest, 
                                                  numBlk = numBlk, numRep = numRep,  
                                                  genLayout = genLayout, 
                                                  numRowPerBlk = numRowPerBlk, numColPerBlk = numColPerBlk,
                                                  numRowPerRep = numRowPerRep, numColPerRep = numColPerRep,
                                                  numFieldRow = numFieldRow, numFieldCol = numFieldCol,
                                                  serpentine = serpentine, topToBottom = topToBottom, 
                                                  checkTrmtList = checkTrmtList, testTrmtList = testTrmtList)
                
                fieldbook <- rbind(fieldbook, data.frame(Trial = i, tmpResult$fieldbook))
                
                if (!genLayout) {
                        trmtLayout <- rbind(trmtLayout, data.frame(Trial = i, tmpResult$plan[[1]]))
                        blkLayout <- rbind(blkLayout, data.frame(Trial = i, tmpResult$plan[[2]]))
                } else {
                        if (i == 1) {
                                trmtLayout <- list()
                                blkLayout <- tmpResult$plan[[4]]
                                plotNumLayout <- tmpResult$plan[[2]]
                                repLayout <- tmpResult$plan[[3]] }
                        trmtLayout[[i]] <- tmpResult$plan[[1]]
                }
        }
        
        if (!genLayout) {
                
                fieldbook <- fieldbook[,c("Trial", "REP", "Block", "ID", "PlotNumber")]
                names(fieldbook) <- c("Trial", "Rep", "Block", "Entry", "PlotNumber")
                entryList <- data.frame(Entry = c(testTrmtList,checkTrmtList), EntryType = c(rep("Test",numTest),rep("Check",numCheck)))
                fieldbook <- merge(fieldbook,entryList)
                fieldbook <- fieldbook[,c("Trial","Rep","Block","Entry","EntryType","PlotNumber")]
                fieldbook <- fieldbook[order(fieldbook$Trial, fieldbook$PlotNumber), ]
                rownames(fieldbook) <- 1:nrow(fieldbook)
                
                if(!treatList){
                  fieldbook <- fieldbook[,c("Trial","Rep","Block","Entry","PlotNumber")]
                } else {
                  names(fieldbook) <- c("Trial", "Rep", "Block", trmtName, paste(trmtName,"Type",sep=""), "PlotNumber")
                }
                
                exptLayout <- list()
                exptLayout[[1]] <- trmtLayout
                exptLayout[[2]] <- blkLayout
                names(exptLayout) <- c("TrmtLayout", "BlockLayout")
                
        } else {
                
                exptLayout <- list()
                exptLayout[[1]] <- trmtLayout
                names(exptLayout[[1]]) <- paste("Trial", 1:length(trmtLayout), sep = "")
                exptLayout[[2]] <- plotNumLayout
                exptLayout[[3]] <- repLayout
                exptLayout[[4]] <- blkLayout
                names(exptLayout) <- c("TrmtLayout", "PlotNumLayout", "RepLayout", "BlockLayout")
                
                tmpResult2 <- list(fieldbook = fieldbook, plan = exptLayout)
                
                tmpResult2$fieldbook <- tmpResult2$fieldbook[,c("Trial", "REP", "Block", "ID", "PlotNumber", "ROW","RANGE")]
                names(tmpResult2$fieldbook) <- c("Trial", "Rep", "Block", "Entry", "PlotNumber", "FieldRow", "FieldColumn")
                entryList <- data.frame(Entry = c(testTrmtList,checkTrmtList), EntryType = c(rep("Test",numTest),rep("Check",numCheck)))
                tmpResult2$fieldbook <- merge(tmpResult2$fieldbook,entryList)
                tmpResult2$fieldbook <- tmpResult2$fieldbook[,c("Trial", "Rep", "Block", "Entry", "EntryType", "PlotNumber", "FieldRow", "FieldColumn")]

                if(!treatList){
                  tmpResult2$fieldbook <- tmpResult2$fieldbook[,c("Trial", "Rep", "Block", "Entry", "PlotNumber", "FieldRow", "FieldColumn")]
                } else {
                  names(tmpResult2$fieldbook) <- c("Trial", "Rep", "Block", trmtName, paste(trmtName,"Type",sep=""), "PlotNumber", "FieldRow", "FieldColumn")
                }
                
                tmpResult2$fieldbook <- tmpResult2$fieldbook[order(tmpResult2$fieldbook$Trial, tmpResult2$fieldbook$PlotNumber),]
                rownames(tmpResult2$fieldbook) <- 1:nrow(tmpResult2$fieldbook)
                
                for (i in 1:length(tmpResult2$plan$TrmtLayout)) {
                        dimnames(tmpResult2$plan$TrmtLayout[[i]]) <- list(paste("FieldRow",1:nrow(tmpResult2$plan$TrmtLayout[[i]]), sep = ""),
                                                                          paste("FieldColumn",1:ncol(tmpResult2$plan$TrmtLayout[[i]]), sep = "")) 
                }
                for (i in 2:length(tmpResult2$plan)) {
                        dimnames(tmpResult2$plan[[i]]) <- list(paste("FieldRow",1:nrow(tmpResult2$plan[[i]]), sep = ""),
                                                               paste("FieldColumn",1:ncol(tmpResult2$plan[[i]]), sep = "")) 
                }
        } ## end stmt -- if (!genLayout) - else stmt
        
        if (display) {
                cat(toupper("Design Properties:"), "\n", sep = "")
                cat("\t", "Incomplete Block Design", "\n", sep = "")
                cat("\t", "Augmented Alpha Lattice Design", "\n\n", sep = "")
                cat(toupper("Design Parameters"), "\n", sep = "")
                cat("\t", "Number of Trials = ", numTrial, "\n", sep = "")
                cat("\t", "Number of Replicated Treatments = ", numCheck, "\n", sep = "")
                cat("\t", "Levels of Replicated Treatments = ", sep = "")
                if (numCheck <= 5) { cat(paste(checkTrmtList, collapse = ", ", sep = ""), "\n", sep = "") 
                } else {
                        cat(paste(checkTrmtList[1:3], collapse = ", ", sep = ""), sep = "") 
                        cat(", ..., ", checkTrmtList[numCheck], "\n", sep = "") 
                }
                cat("\t", "Number of Replicates = ", numRep, "\n", sep = "")
                cat("\t", "Number of UnReplicated Treatments = ", numTest, "\n", sep = "")
                cat("\t", "Levels of UnReplicated Treatments = ", sep = "")
                if (numTest <= 5) { cat(paste(testTrmtList, collapse = ", ", sep = ""), "\n", sep = "") 
                } else {
                        cat(paste(testTrmtList[1:3], collapse = ", ", sep = ""), sep = "") 
                        cat(", ..., ", testTrmtList[numTest], "\n", sep = "") 
                }
                cat("\t", "Number of Plots per Block = ", blksize,"\n", sep = "")
                cat("\t", "Number of Blocks per Replicate = ", numPlotPerRep/blksize,"\n\n", sep = "")
                
                if (genLayout) {
                        cat("\t","Number of Field Rows = ", numFieldRow, "\n",sep = "")
                        cat("\t","Number of Field Columns = ", numFieldCol, "\n\n",sep = "")
                }
                
        }
        
        if (blksize > 30) { cat(green("WARNING: The block size is too large. You might want to consider a smaller block size.","\n\n"))}
        
        if (!genLayout) { 
                return(invisible(list(fieldbook = fieldbook, plan = exptLayout)))  
        } else { 
                return(invisible(list(fieldbook = tmpResult2$fieldbook, plan = tmpResult2$plan))) }
        
}  ## end -- designAugmentedAL function
