# -------------------------------------------------------------------------------------
# Name             : designAlphaLatticeDiGGer
# Description      : Generate randomization for Alpha Lattice Design using DiGGer package
# R Version        : 4.1.0
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 4.1.3
# Modified by      : Justine B. Bonifacio
# Date Modified    : 2021.09.09
# Changes made     : use function from DiGGer package only
# -------------------------------------------------------------------------------------
#' @title Alpha Lattice Design using DiGGer package
#' @aliases designAlphaLatticeDiGGer
#  @aliases designAlphaLatticeDiGGer.default
#'  
#' @description Generate randomization and layout for Alpha Lattice Design using 
#'              DiGGer package.
#'
#' @param generate list of entries to be randomized
#' @param numBlk number of blocks
#' @param numRep number of replicates
#' @param numTrial number of trials
#' @param genLayout logical, whether a layout of the design will be generated
#' @param numRowPerBlk number of rows per block
#' @param numRowPerRep number of rows per replicate
#' @param numFieldRow number of field rows
#' @param serpentine a logical variable indicating whether plot number will be arranged as serpentine order
#' @param topToBottom a logical variable indicating whether plot number will be written from to to bottom
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
#' ## Generate randomization for an experiment with 24 treatment levels in Alpha Lattice replicated 4 times.
#' ## Each replicate will have 4 blocks.
#' alpha1 <- designAlphaLatticeDiGGer(generate = list(Entry = 24), numBlk = 4, numRep = 4, numTrial = 1, 
#'                                    genLayout = FALSE)
#' 
#' ## Generate randomization and layout for an experiment with 24 levels in Alpha Lattice replicated 4 times.
#' ## Each replicate will have 4 blocks. The experiment will be arrange in a 8 x 12 field.
#' ## Each replicate will be arrange in a 4 x 6, while each block will be arrange in 2 x 3.
#' alpha2 <- designAlphaLatticeDiGGer(generate = list(Entry = 24), numBlk = 4, numRep = 4, numTrial = 1, 
#'                                    genLayout = TRUE, numRowPerBlk = 2, numRowPerRep = 4, numFieldRow = 8)
#' 
# -------------------------------------------------------------------------------------

#library(crayon)

designAlphaLatticeDiGGer <- function(generate, numBlk = 2, numRep = 2, numTrial = 1, 
                                     genLayout = FALSE, numRowPerBlk = 1, 
                                     numRowPerRep = 1, numFieldRow = 1,
                                     serpentine = FALSE, topToBottom = TRUE,
                                     display = TRUE){ # UseMethod("designAlphaLatticeDiGGer")}

# designAlphaLatticeDiGGer.default <- function(generate, numBlk = 2, numRep = 2, numTrial = 1, 
#                                              genLayout = FALSE, numRowPerBlk = 1, numRowPerRep = 1, numFieldRow = 1,
#                                              serpentine = FALSE, display = TRUE) {
  
  # --- check inputs --- #
  
  if (numRep < 2) { 
    stop("ERROR: The number of replicates should be greater than or equal to 2.")} # check number of replicates, should be 
  
  if (numBlk == 1) { 
    stop("ERROR: The number of blocks per replicate should be greater than 1.") }
  
  if (length(generate[[1]]) == 1) { 
    tempComb <- FactorList(generate) 
  } else { 
    tempComb <- generate }
  
  flag <- primeNumber(length(tempComb[[1]]))
  if(flag == 1) {
    stop("ERROR: The number of treatments should not be a prime number.")}
  
  if(length(tempComb[[1]]) < 6){
    stop("ERROR: The number of treatments should be greater than or equal to 6.")}
  
  blksize <- length(tempComb[[1]])/numBlk
  if (blksize == 1) { 
    stop("ERROR: The block size should be greater than 1.") }
  
  if (!genLayout) {
    numFieldRow <- length(tempComb[[1]])
    numRowPerRep <- length(tempComb[[1]])
    numRowPerBlk <- blksize
  }
  
  # determine the total number of experimental units
  if ((numRep * length(tempComb[[1]])) > 1500) { 
    stop("ERROR: The maximum number of experimental units that can be generated is 1500.") }
  
  # check if the number of treatment is divisible by the block size
  if (length(tempComb[[1]])%%blksize != 0) { 
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
  if (length(tempComb[[1]])%%numRowPerRep != 0) { 
    stop("ERROR: The number of treatments should be divisible by the number of rows in each replicate.") }
  
  numColPerRep <- length(tempComb[[1]])/numRowPerRep           # determine the number of columns per replicate
  
  if (numFieldRow%%numRowPerRep != 0) { 
    stop("ERROR: The total number of plots should be divisible by the number of field rows.") }
  
  numRepRow <- numFieldRow/numRowPerRep                     # determine the number of rep along the length of the field layout
  
  if (!(numRepRow %in% allFactors(numRep))) { 
    stop("ERROR: The quotient of the number of field rows and number of rows in each replicate should be a factor of the number of the replicates.") }
  
  if((length(tempComb[[1]])*numRep)%%numFieldRow != 0) { 
    stop("ERROR: The total number of plots should be divisible by the number of field rows.") } 
  
  numFieldCol <- (length(tempComb[[1]])*numRep)/numFieldRow    # determine the number of field column in the experiment
  
  fieldbook <- NULL
  trmtLayout <- NULL
  blkLayout <- NULL
  tmpResult2 <- NULL
  
  if (numRep * length(tempComb[[1]]) != numFieldRow * numFieldCol) { stop("Total of plots cannot be accomodated in the field experiment.") }  
  
  for (i in (1:numTrial)) {
    tmpResult <- randomizeIBD(numTrmt = length(tempComb[[1]]), 
                              numFieldRow, numFieldCol, 
                              numRep, numRowPerRep, numColPerRep,
                              numBlk, blksize, numRowPerBlk, numColPerBlk,
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
    fieldbook$Block <- rep(rep(c(1:numBlk), each = blksize), times = numTrial)
    fieldbook$PlotNumber <- as.numeric(paste(fieldbook$REP, paste(rep(0,nchar(length(tempComb[[1]]))), collapse=""), sep = "")) + 1:length(tempComb[[1]]) 
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
    tmpResult2 <- layoutIBD(fieldbook, trmtLayout,
                            numFieldRow, numFieldCol, 
                            numRowPerRep, numColPerRep, 
                            numBlk, numRowPerBlk, numColPerBlk,
                            serpentine = serpentine,
                            topToBottom = topToBottom)
    
    tmpResult2$fieldbook <- tmpResult2$fieldbook[,c("Trial", "REP", "Block", "ID", "PlotNumber", "ROW","RANGE")]
    names(tmpResult2$fieldbook) <- c("Trial", "Rep", "Block", names(tempComb)[1],"PlotNumber", "FieldRow", "FieldColumn")
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
    cat("\t","Alpha Lattice Design","\n\n",sep = "") 
    cat(toupper("Design Parameters:"),"\n",sep = "")
    cat("\t","Number of Trials = ", numTrial, "\n",sep = "")
    cat("\t","Number of Treatments = ", length(tempComb[[1]]), "\n",sep = "")
    cat("\t","Number of Replicates = ", numRep, "\n",sep = "")
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
  
} ## end -- designAlphaLatticeDiGGer function
