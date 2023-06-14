# -------------------------------------------------------------------------------------
# Name             : designAlphaLatticeWithDiagCheck
# Description      : Generate randomization for Alpha Lattice Design with Diagonal Checks
# R Version        : 4.1.2
# -------------------------------------------------------------------------------------
# Author           : Justine B. Dayrit
# Author Email     : j.bonifacio@irri.org
# Script Version   : 1.0.0
# Modified by      : 
# Date Modified    : 
# Changes made     : 
# -------------------------------------------------------------------------------------
#' @title Alpha Lattice Design with Diagonal Checks
#' @aliases designAlphaLatticeWithDiagCheck
#  @aliases designAlphaLatticeWithDiagCheck.default
#'  
#' @description Generate randomization and layout for Alpha Lattice Design with Diagonal Checks.
#'
#' @param numTrmt number of entries to be randomized
#' @param numSpatCheck number of spatial checks
#' @param numBlk number of blocks
#' @param numRep number of replicates
#' @param numTrial number of trials (randomization set-ups with the complete set of entries)
#' @param numRowPerBlk number of rows per block
#' @param numRowPerRep number of rows per replicate
#' @param numFieldRow number of field rows
#' @param treatName NULL or vector of names associated with the treatment
#' @param listName column name associated with the list provided (treatName)
#' @param topToBottom a logical variable indicating whether plot number will be written from to to bottom
#' @param serpentine a logical variable indicating whether plot number will be arranged as serpentine order
#' @param display a logical variable indicating whether randomization parameters will be displayed
#' 
#' @details 
#' Parameters \code{numFieldRow}, \code{numRowPerRep}  and \code{numRowPerBlk} 
#' must be specified.  Values of \code{numRowPerBlk} should be a factor of \code{numRowPerRep} while values of 
#' \code{numRowPerRep} should be a factor of \code{numFieldRow}.
#' 
#' @return A list containing the following components:
#' \item{fieldbook}{a data frame}
#' \item{plan}{a list containing the following components:}
#' \item{TrmtLayout}{a list whose length is equal to the 
#'      number of trials containing the treatment layout for each trial} 
#' \item{PlotNumLayout}{a matrix containing the plot number layout of the experiment} 
#' \item{RepLayout}{a matrix containing the replication layout of the experiment} 
#' \item{BlockLayout}{a matrix containing the block layout of the experiment}
#' 
#' @examples
#' ## Generate randomization and layout for an experiment with 40 entries replicated 2 times.
#' ## Each replicate will have 4 blocks. The experiment will be arrange in a 4 x 20 field.
#' ## Each replicate will be arrange in a 2 x 20, while each block will be arrange in 1 x 10.
#' ## Inserting 1 spatial check in each block diagonally.
#' alphadc1 <- designAlphaLatticeWithDiagCheck (numTrmt = 24, numSpatCheck = 1,
#'                                             numBlk = 4, numRep = 2, numTrial = 1,
#'                                             numRowPerBlk = 1, numRowPerRep = 2,
#'                                             numFieldRow = 4, serpentine = FALSE, topToBottom = FALSE,
#'                                             display = TRUE)
#' 
# -------------------------------------------------------------------------------------

#library(crayon)

designAlphaLatticeWithDiagCheck <- function(numTrmt, numSpatCheck = 1, numBlk = 2, numRep = 2, 
                                            numTrial = 1, numRowPerBlk = 1, numRowPerRep = 1,
                                            numFieldRow = 1, serpentine = FALSE, topToBottom = TRUE,
                                            treatName = NULL, listName = NULL,
                                            display = TRUE) { #UseMethod("designAlphaLatticeWithDiagCheck") }
  
  
  #designAlphaLatticeWithDiagCheck.default <- function(numTrmt, numSpatCheck = 1, numBlk = 2, numRep = 2, 
  #                                                    numTrial = 1, numRowPerBlk = 1, numRowPerRep = 1,
  #                                                    numFieldRow = 1, serpentine = FALSE, topToBottom = TRUE,
  #                                                    display = TRUE) { 
  
  # --- check inputs --- #
  
  if (numRep < 2) { 
    stop("ERROR: The number of replicates should be greater than or equal to 2.")}
  
  if (numBlk == 1) { 
    stop("ERROR: The number of blocks per replicate should be greater than 1.") }
  
  if (numSpatCheck < 1) { 
    stop("ERROR: There should be at least 1 spatial check.")}
  
  if (numSpatCheck > (numBlk*numRep)) { 
    stop("ERROR: The number of spatial checks should be less than the total number of blocks.") }
  
  treatList <- TRUE
  # default value when treatName is NULL
  if(is.null(treatName)){
    treatName <- c(paste0("Entry",1:numTrmt),paste0("SpatialCheck",1:numSpatCheck))
    treatList <- FALSE
  }
  
  flag <- primeNumber(numTrmt)
  if(flag == 1) {
    stop("ERROR: The number of treatments should not be a prime number.")}
  
  if(numTrmt < 6){
    stop("ERROR: The number of treatments should be greater than or equal to 6.")}
  
  blksize <- numTrmt/numBlk
  if (blksize == 1) { 
    stop("ERROR: The block size should be greater than 1.") }
  
  # check if the number of treatment is divisible by the block size
  if (numTrmt%%blksize != 0) { 
    stop("ERROR: The number of treatments should be divisible by the number of plots 
         per block.") }
  
  if (blksize%%numRowPerBlk != 0) { 
    stop("ERROR: The number of plots per block should be divisible by the number of 
         rows per block.") }
  
  # determine the number of columns per block within a replicate
  numColPerBlk <- blksize/numRowPerBlk                         
  
  # check if number of rows per replicate is divisible by the number of rows per block
  if (numRowPerRep%%numRowPerBlk != 0) { 
    stop("ERROR: The number of rows per replicate should be divisible by the number of 
         rows per block.") }
  
  # check if the quotient of number of rows per replicate and number of rows per block 
  # is a factor to number of blocks per replicate
  if (!((numRowPerRep/numRowPerBlk) %in% allFactors(numBlk))) { 
    stop("ERROR: The quotient of the number of rows in each replicate and number of rows 
         in each block should be a factor of the number of blocks per replicate.") }
  
  # check if the treatment is divisible by the number of rows per rep
  if (numTrmt%%numRowPerRep != 0) { 
    stop("ERROR: The number of treatments should be divisible by the number of rows in 
         each replicate.") }
  
  # determine the number of columns per replicate
  numColPerRep <- numTrmt/numRowPerRep           
  
  if (numFieldRow%%numRowPerRep != 0) { 
    stop("ERROR: The total number of plots should be divisible by the number of 
         field rows.") }
  
  if (numFieldRow == 1) { 
    stop("ERROR: The number of field rows should be greater than 1.") }
  
  # determine the number of rep along the width of the field layout
  numRepRow <- numFieldRow/numRowPerRep                     
  
  if (!(numRepRow %in% allFactors(numRep))) { 
    stop("ERROR: The quotient of the number of field rows and number of rows in 
         each replicate should be a factor of the number of the replicates.") }
  
  if((numTrmt*numRep)%%numFieldRow != 0) { 
    stop("ERROR: The total number of plots should be divisible by the number of 
         field rows.") } 
  
  # determine the number of field column in the experiment
  numFieldCol <- (numTrmt*numRep)/numFieldRow    
  
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

  if (numRep * numTrmt != numFieldRow * numFieldCol) { 
    stop("Total of plots cannot be accomodated in the field experiment.") }
  
  for (i in (1:numTrial)) {
    
    if((numRep == 2 & blksize <= numBlk) | 
       (numRep == 3 & numBlk%%2 != 0 & blksize <= numBlk) |
       (numRep == 3 & numBlk%%2 == 0 & blksize <= numBlk-1) |
       (numRep == 4 & numBlk%%2 != 0 & numBlk%%3 != 0 & blksize <= numBlk)){
      
      tmpResult <- randomizeAlphaLatticeWithDiagCheck(numTrmt = numTrmt, numSpatCheck = numSpatCheck,
                                                      numFieldRow = numFieldRow, numFieldCol = numFieldCol,
                                                      numRep = numRep, numRowPerRep = numRowPerRep, numColPerRep = numColPerRep,
                                                      numBlk = numBlk, blksize = blksize, numRowPerBlk = numRowPerBlk, numColPerBlk = numColPerBlk,
                                                      topToBottom = topToBottom, serpentine = serpentine, trmtList = treatName)
      
      tmpResult$designparam <- "Agricolae Alpha-Latice (0,1)"
      
    } else{
      tmpResult <- randomizeAlphaLatticeEBSWithDiagCheck(numTrmt = numTrmt, numSpatCheck = numSpatCheck,
                                                         numFieldRow = numFieldRow, numFieldCol = numFieldCol,
                                                         numRep = numRep, numRowPerRep = numRowPerRep, numColPerRep = numColPerRep,
                                                         numBlk = numBlk, blksize = blksize, numRowPerBlk = numRowPerBlk, numColPerBlk = numColPerBlk,
                                                         topToBottom = topToBottom, serpentine = serpentine, trmtList = treatName)
      
      tmpResult$designparam <- paste("ebsRtools",tmpResult$designparam,sep=" ")
    }

    fieldbook <- rbind(fieldbook, data.frame(Trial = i, tmpResult$fieldbook))
    
    if (i == 1) {
      trmtLayout <- list()
      blkLayout <- tmpResult$plan[[4]]
      plotNumLayout <- tmpResult$plan[[2]]
      repLayout <- tmpResult$plan[[3]] }
    trmtLayout[[i]] <- tmpResult$plan[[1]]
  }
  
  exptLayout <- list()
  exptLayout[[1]] <- trmtLayout
  names(exptLayout[[1]]) <- paste("Trial", 1:length(trmtLayout), sep = "")
  exptLayout[[2]] <- plotNumLayout
  exptLayout[[3]] <- repLayout
  exptLayout[[4]] <- blkLayout
  names(exptLayout) <- c("TrmtLayout", "PlotNumLayout", "RepLayout", "BlockLayout")
  
  tmpResult2 <- list(fieldbook = fieldbook, plan = exptLayout)
  
  if(treatList){
    tmpResult2$fieldbook <- mutate(tmpResult2$fieldbook, TYPE = if_else(as.numeric(tmpResult2$fieldbook$ENTRY) <= numTrmt, "Entry", "Spatial_Check"))
    tmpResult2$fieldbook <- tmpResult2$fieldbook[,c("Trial", "REP", "Block", "ID", "TYPE", "PlotNumber", "ROW","RANGE")]
    if(!is.null(listName)){
      names(tmpResult2$fieldbook) <- c("Trial", "Rep", "Block", listName, paste0(listName,"Type"), "PlotNumber", "FieldRow", "FieldColumn")
    } else{
      names(tmpResult2$fieldbook) <- c("Trial", "Rep", "Block", "Entry", "EntryType", "PlotNumber", "FieldRow", "FieldColumn")
    }
  } else{
    tmpResult2$fieldbook <- tmpResult2$fieldbook[,c("Trial", "REP", "Block", "ID", "PlotNumber", "ROW","RANGE")]
    if(!is.null(listName)){
      names(tmpResult2$fieldbook) <- c("Trial", "Rep", "Block", listName, "PlotNumber", "FieldRow", "FieldColumn")
    } else{
      names(tmpResult2$fieldbook) <- c("Trial", "Rep", "Block", "Entry", "PlotNumber", "FieldRow", "FieldColumn")
    }
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
  
  if (display) {
    cat(toupper("Design Properties:"),"\n",sep = "")
    cat("\t","Incomplete Block Design","\n",sep = "") 
    cat("\t","Alpha Lattice Design with Diagonal Checks","\n\n",sep = "") 
    cat(toupper("Design Parameters:"),"\n",sep = "")
    cat("\t","Number of Trials = ", numTrial, "\n",sep = "")
    cat("\t","Number of Treatments = ", numTrmt, "\n",sep = "")
    cat("\t","Number of Spatial Checks = ", numSpatCheck, "\n",sep = "")
    cat("\t","Number of Replicates = ", numRep, "\n",sep = "")
    cat("\t","Number of Plots per Block = ", blksize+1, "\n",sep = "")
    cat("\t","Number of Blocks per Replicate = ", numBlk, "\n\n",sep = "")
    
    if(numRowPerBlk == 1){
      numBlkPerRow <- numFieldCol/numColPerBlk
      numFieldCol <- numFieldCol + numBlkPerRow
    } else if(numRowPerBlk == blksize){
      numBlkPerCol <- numFieldRow/numRowPerBlk
      numFieldRow <- numFieldRow + numBlkPerCol
    }
    
    cat("\t","Number of Field Rows = ", numFieldRow, "\n",sep = "")
    cat("\t","Number of Field Columns = ", numFieldCol, "\n\n",sep = "")
    
    cat("\t",tmpResult$designparam, "\n",sep = "")
    
  }
  
  if (blksize > 30) { cat(green("WARNING: The block size is too large for an Alpha Lattice. You might want to consider a smaller block size.","\n\n"))}
  
  return(invisible(list(fieldbook = tmpResult2$fieldbook, plan = tmpResult2$plan, concurrence = tmpResult$designparam)))

} ## end -- designAlphaLatticeWithDiagCheck function
