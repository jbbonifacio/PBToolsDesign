# -------------------------------------------------------------------------------------
# Name             : designAlphaLattice
# Description      : Generate randomization for Alpha Lattice Design 
#                    using agricolae and ebsRTools package
# R Version        : 4.1.0
# -------------------------------------------------------------------------------------
# Author           : Justine B. Bonifacio 
# Author Email     : j.bonifacio@irri.org
# Script Version   : 1.0.1
# Modified by      : 
# Date Modified    : June 9,2023
# Changes made     : added k!=2 condition for agricolae function
# -------------------------------------------------------------------------------------
#' @title Alpha Lattice Design using agricolae and ebsRTools package
#' @aliases designAlphaLattice
#  @aliases designAlphaLattice.default
#'  
#' @description Generate randomization and layout for Alpha Lattice Design using 
#'              agricolae and ebsRTools package.
#'
#' @param generate list of entries to be randomized
#' @param numBlk number of blocks
#' @param numRep number of replicates
#' @param numTrial number of trials
#' @param genLayout logical, whether a layout of the design will be generated
#' @param numRowPerBlk number of rows per block
#' @param numRowPerRep number of rows per replicate
#' @param numFieldRow number of field rows
#' @param topToBottom logical, whether plot number will be written from to to bottom
#' @param serpentine logical, whether plot number will be arranged as serpentine order
#' @param display logical, whether randomization parameters will be displayed
#' 
#' @details 
#' If \code{genLayout} is \code{TRUE}, then parameters \code{numFieldRow}, 
#' \code{numRowPerRep}  and \code{numRowPerBlk} must be specified.  
#' Values of \code{numRowPerBlk} should be a factor of \code{numRowPerRep} 
#' while values of \code{numRowPerRep} should be a factor of \code{numFieldRow}.
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
#' ## Generate randomization for an experiment with 24 treatment levels in 
#' ## Alpha Lattice replicated 4 times. Each replicate will have 4 blocks.
#' alpha1 <- designAlphaLattice(generate = list(Entry = 24), numBlk = 4, numRep = 4,  
#'                              numTrial = 1, genLayout = FALSE)
#' 
#' ## Generate randomization and layout for an experiment with 24 levels in 
#' ## Alpha Lattice replicated 4 times. Each replicate will have 4 blocks. 
#' ## The experiment will be arrange in a 8 x 12 field.
#' ## Each replicate will be arrange in a 4 x 6, while each block will be arrange in 2 x 3.
#' alpha2 <- designAlphaLattice(generate = list(Entry = 24), numBlk = 4, numRep = 4,  
#'                              numTrial = 1, genLayout = TRUE, numRowPerBlk = 2, 
#'                              numRowPerRep = 4, numFieldRow = 8)
#' 
# -------------------------------------------------------------------------------------

# library(crayon)

designAlphaLattice <- function(generate, numBlk = 2, numRep = 2, numTrial = 1,
                               genLayout = FALSE, numRowPerBlk = 1, numRowPerRep = 1,
                               numFieldRow = 1, serpentine = FALSE, topToBottom = TRUE,
                               display = TRUE) { #UseMethod("designAlphaLattice") }
  
# designAlphaLattice.default <- function(generate, numBlk = 2, numRep = 2, numTrial = 1, 
#                                        genLayout = FALSE, numRowPerBlk = 1, numRowPerRep = 1,
#                                        numFieldRow = 1, serpentine = FALSE,
#                                        topToBottom = TRUE, display = TRUE) {

    
  # --- check inputs --- #
  
  if (numRep < 2) {
    stop("ERROR: The number of replicates should be greater than or equal to 2.")}
  
  if (numBlk == 1) {
    stop("ERROR: The number of blocks per replicate should be greater than 1.")}
  
  if (length(generate[[1]]) == 1) { 
    tempComb <- FactorList(generate) 
  } else { 
    tempComb <- generate}
  
  flag <- primeNumber(length(tempComb[[1]]))
  if(flag == 1){
    stop("ERROR: The number of treatments should not be a prime number.")}
  
  if(length(tempComb[[1]]) < 6) {
    stop("ERROR: The number of treatments should be greater than or equal to 6.")}
  
  blksize <- length(tempComb[[1]])/numBlk
  if (blksize == 1) { 
    stop("ERROR: The block size should be greater than 1.")}
  
  if (!genLayout) {
    numFieldRow <- length(tempComb[[1]])
    numRowPerRep <- length(tempComb[[1]])
    numRowPerBlk <- blksize
  }
  
  # check if the number of treatment is divisible by the block size
  if (length(tempComb[[1]])%%blksize != 0) { 
    stop("ERROR: The number of treatments should be divisible by the number of plots 
         per block.")}

  if (blksize%%numRowPerBlk != 0) { 
    stop("ERROR: The number of plots per block should be divisible by the number of 
         rows per block.")}
  
  # determine the number of columns per block within a replicate
  numColPerBlk <- blksize/numRowPerBlk                         
  
  # check if number of rows per replicate is divisible by the number of rows per block
  if (numRowPerRep%%numRowPerBlk != 0) { 
    stop("ERROR: The number of rows per replicate should be divisible by the number of 
         rows per block.")}
  
  # check if the quotient of number of rows per replicate and number of rows per block 
  # is a factor to number of blocks per replicate
  if (!((numRowPerRep/numRowPerBlk) %in% allFactors(numBlk))) { 
    stop("ERROR: The quotient of the number of rows in each replicate and number of rows 
         in each block should be a factor of the number of blocks per replicate.")}
  
  # check if the treatment is divisible by the number of rows per rep
  if (length(tempComb[[1]])%%numRowPerRep != 0) { 
    stop("ERROR: The number of treatments should be divisible by the number of rows in 
         each replicate.")}
  
  # determine the number of columns per replicate
  numColPerRep <- length(tempComb[[1]])/numRowPerRep           
  
  if (numFieldRow%%numRowPerRep != 0) { 
    stop("ERROR: The total number of plots should be divisible by the number of 
         field rows.")}
  
  # determine the number of rep along the width of the field layout
  numRepRow <- numFieldRow/numRowPerRep                     
  
  if (!(numRepRow %in% allFactors(numRep))) { 
    stop("ERROR: The quotient of the number of field rows and number of rows in 
         each replicate should be a factor of the number of the replicates.")}
  
  if((length(tempComb[[1]])*numRep)%%numFieldRow != 0) { 
    stop("ERROR: The total number of plots should be divisible by the number of 
         field rows.")} 
  
  # determine the number of field column in the experiment
  numFieldCol <- (length(tempComb[[1]])*numRep)/numFieldRow    
  
  # determine the number of replicates along the length of the field layout
  numRepCol <- numFieldCol/numColPerRep                     
  # determine the number of blocks along the width of each replicate
  numBlkRow <- numRowPerRep/numRowPerBlk 
  # determine the number of blocks along the length of each replicate
  numBlkCol <- numColPerRep/numColPerBlk                       
  
  fieldbook <- NULL
  trmtLayout <- NULL
  blkLayout <- NULL
  plotNumLayout <- NULL
  repLayout <- NULL
  tmpResult2 <- NULL
  
  if (numRep * length(tempComb[[1]]) != numFieldRow * numFieldCol) { 
    stop("Total of plots cannot be accomodated in the field experiment.")}
  
  for (i in (1:numTrial)) {
    
    if((numRep == 2 & blksize <= numBlk & blksize != 2) | 
       (numRep == 3 & numBlk%%2 != 0 & blksize <= numBlk) |
       (numRep == 3 & numBlk%%2 == 0 & blksize <= numBlk-1 & blksize != 2) |
       (numRep == 4 & numBlk%%2 != 0 & numBlk%%3 != 0 & blksize <= numBlk)) {
      
      tmpResult <- randomizeAlphaLattice(numTrmt = length(tempComb[[1]]), 
                                         numFieldRow = numFieldRow, numFieldCol = numFieldCol,
                                         numRep = numRep, numRowPerRep = numRowPerRep, numColPerRep = numColPerRep,
                                         numBlk = numBlk, blksize = blksize, numRowPerBlk = numRowPerBlk, numColPerBlk = numColPerBlk,
                                         numBlkRow = numBlkRow, numBlkCol = numBlkCol, numRepRow = numRepRow, numRepCol = numRepCol,
                                         topToBottom = topToBottom, serpentine = serpentine, trmtList = tempComb[[1]],
                                         genLayout = genLayout)
      
      tmpResult$designparam <- "Agricolae Alpha-Latice (0,1)"
    
    } else {
      tmpResult <- randomizeAlphaLatticeEBS(numTrmt = length(tempComb[[1]]), 
                                            numFieldRow = numFieldRow, numFieldCol = numFieldCol,
                                            numRep = numRep, numRowPerRep = numRowPerRep, numColPerRep = numColPerRep,
                                            numBlk = numBlk, blksize = blksize, numRowPerBlk = numRowPerBlk, numColPerBlk = numColPerBlk,
                                            numBlkRow = numBlkRow, numBlkCol = numBlkCol, numRepRow = numRepRow, numRepCol = numRepCol,
                                            topToBottom = topToBottom, serpentine = serpentine, trmtList = tempComb[[1]],
                                            genLayout = genLayout)
      
      tmpResult$designparam <- paste("ebsRtools",tmpResult$designparam,sep=" ")
    }
    
    
    fieldbook <- rbind(fieldbook, data.frame(Trial = i, tmpResult$fieldbook))
    
    if (!genLayout) {
      trmtLayout <- rbind(trmtLayout, data.frame(Trial = i, tmpResult$plan[[1]]))
      blkLayout <- rbind(blkLayout, data.frame(Trial = i, tmpResult$plan[[2]]))
    } else {
      if (i == 1) {
        trmtLayout <- list()
        blkLayout <- tmpResult$plan[[4]]
        plotNumLayout <- tmpResult$plan[[2]]
        repLayout <- tmpResult$plan[[3]]}
      trmtLayout[[i]] <- tmpResult$plan[[1]]
    }
  }
  
  if (!genLayout) {
    fieldbook <- fieldbook[,c("Trial", "REP", "Block", "ID", "PlotNumber")]
    names(fieldbook) <- c("Trial", "Rep", "Block", names(tempComb)[1], "PlotNumber")
    
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
    
    cat("\t",tmpResult$designparam, "\n",sep = "")
    
  }
  
  if (blksize > 30) { cat(green("WARNING: The block size is too large for an Alpha Lattice. You might want to consider a smaller block size.","\n\n"))}
  
  if (!genLayout) { 
    return(invisible(list(fieldbook = fieldbook, plan = exptLayout, concurrence = tmpResult$designparam)))  
  } else { 
    return(invisible(list(fieldbook = tmpResult2$fieldbook, plan = tmpResult2$plan, concurrence = tmpResult$designparam))) }
  
} ## end -- designAlphaLattice function
