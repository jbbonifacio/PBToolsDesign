# -------------------------------------------------------------------------------------
# Name             : designAlphaLatticeWithRepCheck
# Description      : Generate randomization for Alpha Lattice Design with Repeated Checks
# R Version        : 4.0.3
# -------------------------------------------------------------------------------------
# Author           : Justine B. Bonifacio
# Author Email     : j.bonifacio@irri.org
# Script Version   : 1.0.2
# Modified by      : Justine B. Bonifacio
# Date Modified    : 2021.10.05
# Changes made     : updated numBlk value checker
# -------------------------------------------------------------------------------------
#' @title Alpha Lattice Design with Repeated Checks
#' @aliases designAlphaLatticeWithRepCheck
#  @aliases designAlphaLatticeWithRepCheck.default
#'  
#' @description Generate randomization and layout for Alpha Lattice Design with Repeated Checks.
#'
#' @param generate list of entries to be randomized
#' @param numRepCheck number of entries repeated within a replicate per group
#' @param numBlk number of blocks
#' @param numRep number of replicates
#' @param numRepOfRepCheck number of times an entry will be repeated within a replicate per group
#' @param numTrial number of trials (randomization set-ups with the complete set of entries)
#' @param genLayout logical, whether a layout of the design will be generated
#' @param numRowPerBlk number of rows per block
#' @param numRowPerRep number of rows per replicate
#' @param numFieldRow number of field rows
#' @param topToBottom a logical variable indicating whether plot number will be written from to to bottom
#' @param serpentine a logical variable indicating whether plot number will be arranged as serpentine order
#' @param display a logical variable indicating whether randomization parameters will be displayed
#' 
#' @details 
#' If \code{genLayout} is \code{TRUE}, then parameters \code{numFieldRow}, \code{numRowPerRep}  and \code{numRowPerBlk} 
#' must be specified.  Values of \code{numRowPerBlk} should be a factor of \code{numRowPerRep} while values of 
#' \code{numRowPerRep} should be a factor of \code{numFieldRow}.
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
#' ## Generate randomization for an experiment with 16 test entries and 2 groups of replicated checks replicated 2 times.
#' ## Group 1: 4 entries replicated 3 times
#' ## Group 2: 3 entry replicated 4 times
#' ## Each replicate will have 4 blocks.
#' 
#' alphawr1 <- designAlphaLatticeWithRepCheck (generate = list(Entry = 23), numRepCheck = c(4,3),
#'                                             numBlk = 4, numRep = 2, numRepOfRepCheck = c(3,4), numTrial = 1,
#'                                             genLayout = FALSE)
#' 
#' ## Generate randomization and layout for an experiment 16 test entries and 2 groups of replicated checks replicated 2 times.
#' ## Group 1: 4 entries replicated 3 times
#' ## Group 2: 3 entry replicated 4 times
#' ## Each replicate will have 4 blocks. The experiment will be arrange in a 20 x 4 field.
#' ## Each replicate will be arrange in a 20 x 2, while each block will be arrange in 10 x 1.
#' alphawr2 <- designAlphaLatticeWithRepCheck (generate = list(Entry = 23), numRepCheck = c(4,3),
#'                                             numBlk = 4, numRep = 2, numRepOfRepCheck = c(3,4), numTrial = 1,
#'                                             genLayout = TRUE, numRowPerBlk = 10, numRowPerRep = 20,
#'                                             numFieldRow = 20, serpentine = FALSE, topToBottom = TRUE,
#'                                             display = TRUE)
#' 
# -------------------------------------------------------------------------------------

#library(crayon)

designAlphaLatticeWithRepCheck <- function(generate, numRepCheck = 1, numBlk = 2, numRep = 2, 
                                           numRepOfRepCheck = 2, numTrial = 1,
                                           genLayout = FALSE, numRowPerBlk = 1, numRowPerRep = 1,
                                           numFieldRow = 1, serpentine = FALSE, topToBottom = TRUE,
                                           display = TRUE) { #UseMethod("designAlphaLatticeWithRepCheck") }


#designAlphaLatticeWithRepCheck.default <- function(generate, numRepCheck = 1, numBlk = 2, numRep = 2, 
#                                                   numRepOfRepCheck = 2, numTrial = 1,
#                                                   genLayout = FALSE, numRowPerBlk = 1, numRowPerRep = 1,
#                                                   numFieldRow = 1, serpentine = FALSE, topToBottom = TRUE,
#                                                   display = TRUE) { 
  
  # --- check inputs --- #
  if (length(generate[[1]]) == 1) { 
    tempComb <- FactorList(generate) 
  } else { 
    tempComb <- generate }
  
  if (sum(numRepCheck) >= length(tempComb[[1]])) { 
    stop("ERROR: The number of Repeated Checks should be less than the number of entries.")}
  
  if (numRep < 2 | any(numRepOfRepCheck < 2)) { 
    stop("ERROR: The number of replicates should be greater than or equal to 2.")} # check number of replicates, should be 
  
  numEntries <- length(tempComb[[1]])
  treatRepPerRep <- rep(c(1,numRepOfRepCheck), c(length(tempComb[[1]])-sum(numRepCheck), numRepCheck))
  treatGroup <- rep(c(1,2), c(length(tempComb[[1]])-sum(numRepCheck), sum(numRepCheck)))
  
  if (numBlk == 1 | any(numRepOfRepCheck > numBlk)) { 
    stop("ERROR: The number of blocks per replicate should be greater than 1 and the number of maximum check repeat.") }
  
  flag <- primeNumber(sum(treatRepPerRep))
  if(flag == 1) {
    stop("ERROR: The number of treatments should not be a prime number.")}
  
  if(sum(treatRepPerRep) < 6){
    stop("ERROR: The number of treatments should be greater than or equal to 6.")}
  
  blksize <- sum(treatRepPerRep)/numBlk
  if (blksize == 1) { 
    stop("ERROR: The block size should be greater than 1.") }
  
  if (!genLayout) {
    numFieldRow <- sum(treatRepPerRep)
    numRowPerRep <- sum(treatRepPerRep)
    numRowPerBlk <- blksize
  }
  
  # determine the total number of experimental units
  if ((numRep * sum(treatRepPerRep)) > 1500) { 
    stop("ERROR: The maximum number of experimental units that can be generated is 1500.") }
  
  # check if the number of treatment is divisible by the number of blocks
  if (sum(treatRepPerRep)%%numBlk != 0) { 
    stop("ERROR: The number of treatments should be divisible by the number of plots per block.") }
  
  if (blksize%%numRowPerBlk != 0) { 
    stop("ERROR: The number of plots per block should be divisible by the number of rows per block.") }
  
  numColPerBlk <- blksize/numRowPerBlk                         # determine the number of columns per block with a replicate
  
  # check if # of rows per replicate is divisible by the # of rows per block
  if (numRowPerRep%%numRowPerBlk != 0) { 
    stop("ERROR: The number of rows per replicate should be divisible by the number of rows per block.") }
  
  # check if the quotient of # of rows per replicate and # of rows per block is a factor to number of blocks per replicate
  if (!((numRowPerRep/numRowPerBlk) %in% allFactors(numBlk))) { 
    stop("ERROR: The quotient of the number of rows in each replicate and number of rows in each block should be a factor of the number of blocks per replicate.") }
  
  # check if the treatment is divisible by the number of rows per Rep
  if (sum(treatRepPerRep)%%numRowPerRep != 0) { 
    stop("ERROR: The number of treatments should be divisible by the number of rows in each replicate.") }
  
  numColPerRep <- sum(treatRepPerRep)/numRowPerRep           # determine the number of columns per replicate
  
  if (numFieldRow%%numRowPerRep != 0) { 
    stop("ERROR: The total number of plots should be divisible by the number of field rows.") }
  
  numRepRow <- numFieldRow/numRowPerRep                     # determine the number of rep along the length of the field layout
  
  if (!(numRepRow %in% allFactors(numRep))) { 
    stop("ERROR: The quotient of the number of field rows and number of rows in each replicate should be a factor of the number of the replicates.") }
  
  if((sum(treatRepPerRep)*numRep)%%numFieldRow != 0) { 
    stop("ERROR: The total number of plots should be divisible by the number of field rows.") } 
  
  numFieldCol <- (sum(treatRepPerRep)*numRep)/numFieldRow    # determine the number of field column in the experiment
  
  numRepCol <- numFieldCol/numColPerRep                     # determine the number of blocks along the length of each replicate
  numBlkRow <- numRowPerRep/numRowPerBlk                       # determine the number of blocks along the length of each replicate
  numBlkCol <- numColPerRep/numColPerBlk                       
  
  fieldbook <- NULL
  trmtLayout <- NULL
  blkLayout <- NULL
  plotNumLayout <- NULL
  repLayout <- NULL
  tmpResult2 <- NULL
  
  if (numRep * sum(treatRepPerRep) != numFieldRow * numFieldCol) { stop("Total of plots cannot be accomodated in the field experiment.") }
  
  for (i in (1:numTrial)) {
    tmpResult <- randomizeAlphaLatticeWithRepCheck(numTrmt = numEntries, numFieldRow, numFieldCol,
                                                   numRowPerRep, numColPerRep, 
                                                   numRowPerBlk, numColPerBlk,
                                                   trmtRepPerRep = treatRepPerRep, 
                                                   trmtGroup = treatGroup,
                                                   trmtList = tempComb[[1]])
    
    fieldbook <- rbind(fieldbook, data.frame(Trial = i, tmpResult$fieldbook))
    
    if (!genLayout) {
      trmtLayout <- rbind(trmtLayout, data.frame(Trial = i, tmpResult$plan))
    } else {
      if (i == 1) trmtLayout <- list()
      trmtLayout[[i]] <- tmpResult$plan
    }
    
  } ## end stmt for (i in (1:trial))
  
  if (!genLayout) {
    fieldbook$Block <- rep(rep(rep(c(1:numBlk), each = blksize), times = numRep), times = numTrial)
    fieldbook$PlotNumber <- as.numeric(paste(fieldbook$REP, paste(rep(0,nchar(sum(treatRepPerRep))), collapse=""), sep = "")) + 1:sum(treatRepPerRep) 
    fieldbook <- fieldbook[,c("Trial", "REP", "Block", "ID", "PlotNumber")]
    fieldbook <- fieldbook[order(fieldbook$Trial, fieldbook$PlotNumber),]
    names(fieldbook) <- c("Trial", "Rep", "Block", names(tempComb)[1], "PlotNumber")
    names(trmtLayout)[2:ncol(trmtLayout)] <- paste("Rep", 1:numRep, sep = "")
    blkLayout <- trmtLayout
    
    for (i in 2:ncol(blkLayout)){
      blkLayout[,i] <- rep(1:numBlk, each = blksize, times = numTrial)}
    
    exptLayout <- list()
    exptLayout[[1]] <- trmtLayout
    exptLayout[[2]] <- blkLayout
    names(exptLayout) <- c("TrmtLayout", "BlockLayout")
    
  } else {
    tmpResult2 <- layoutAlphaLatticeWithRepCheck(fieldbook, trmtLayout,
                                                 numFieldRow, numFieldCol,
                                                 numRowPerRep, numColPerRep,
                                                 numBlk, numRowPerBlk, numColPerBlk,
                                                 serpentine = serpentine,
                                                 topToBottom = topToBottom)  
    
    tmpResult2$fieldbook <- tmpResult2$fieldbook[,c("Trial", "REP", "Block", "ID", "PlotNumber", "ROW","RANGE")]
    names(tmpResult2$fieldbook) <- c("Trial", "Rep", "Block", names(tempComb)[1], "PlotNumber", "FieldRow", "FieldColumn")
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
    cat(toupper("Design Properties:"),"\n",sep = "")
    cat("\t","Incomplete Block Design","\n",sep = "") 
    cat("\t","Alpha Lattice Design with Repeated Checks","\n\n",sep = "") 
    cat(toupper("Design Parameters:"),"\n",sep = "")
    cat("\t","Number of Trials = ", numTrial, "\n",sep = "")
    cat("\t","Number of Treatments = ", numEntries, "\n",sep = "")
    cat("\t","Number of Replicates = ", numRep, "\n",sep = "")
    
    for (i in (1:length(numRepCheck))) {
      cat("\t", "Replicated Check Group ", i, ": ", "\n", sep = "")
      cat("\t", "Number of Replicates for Replicated Check Group ", i, " = ",  numRepOfRepCheck[i],"\n", sep = "")
      cat("\t", "Levels of ", names(numRepCheck)[i]," = ", sep = "")
      
      start <- length(tempComb[[1]]) - sum(numRepCheck[i:length(numRepCheck)]) + 1
      end <- start + numRepCheck[i] - 1
      
      if (numRepCheck[[i]] <= 5) { 
        cat(paste(tempComb[[1]][start:end], collapse = ", ", sep = ""), "\n", sep = "") 
      } else {
        cat(paste(tempComb[[1]][start:(start+2)], collapse = ", ", sep = ""), sep = "") 
        cat(", ..., ", paste(tempComb[[1]][end], sep = ""), "\n", sep = "") 
      }
      
      cat("\n")
    }
    
    cat("\t","Number of Plots per Block = ", blksize, "\n",sep = "")
    cat("\t","Number of Blocks per Replicate = ", numBlk, "\n\n",sep = "")
    
    if (genLayout) {
      cat("\t","Number of Field Rows = ", numFieldRow, "\n",sep = "")
      cat("\t","Number of Field Columns = ", numFieldCol, "\n\n",sep = "")
    }
    
  }
  
  if (blksize > 30) { cat(green("WARNING: The block size is too large for an Alpha Lattice. You might want to consider a smaller block size.","\n\n"))}
  
  if (!genLayout) { 
    return(invisible(list(fieldbook = fieldbook, plan = exptLayout)))  
  } else { 
    return(invisible(list(fieldbook = tmpResult2$fieldbook, plan = tmpResult2$plan))) }
  
} ## end -- designAlphaLatticeWithRepCheck function
