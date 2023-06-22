# -------------------------------------------------------------------------------------
# Name             : designRowColumn
# Description      : Generate randomization for Row-Column Design
# R Version        : 4.0.1
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 4.0.0
# Modified by      :
# Date Modified    :
# Changes made     :
# -------------------------------------------------------------------------------------
#' @title Row-Column Design
#' @aliases designRowColumn
#  @aliases designRowColumn.default
#'  
#' @description Generate randomization and layout for Row-Column Design.
#'
#' @param generate list of entries to be randomized
#' @param numRep number of replicates or blocks
#' @param numTrial number of trials (randomization set-ups with the complete set of entries)
#' @param genLayout logical, whether a layout of the design will be generated or not
#' @param numRowBlk number of rows per replicate or number of row blocks per replicate
#' @param numFieldRow number of field rows
#' @param topToBottom logical, whether plot number will be written from to to bottom
#' @param serpentine logical, whether plot number will be arranged as serpentine order
#' @param display a logical variable indicating whether randomization parameters will be displayed
#' 
#' @details 
#' If \code{genLayout} is \code{TRUE}, then parameters \code{numFieldRow}, \code{numRowBlk}  and \code{numRowPerBlk} 
#' must be specified.  Values of \code{numRowPerBlk} should be a factor of \code{numRowBlk} while values of 
#' \code{numRowBlk} should be a factor of \code{numFieldRow}.
#' 
#' @return A list containing the following components:
#' \item{fieldbook}{a data frame}
#' \item{plan}{a list containing the following components: TrmtLayout, RowBlockLayout and ColumnBlockLayout, 
#' if \code{genLayout} is \code{FALSE} or a list containing the following components 
#' if \code{genLayout} is \code{TRUE}:}
#' \item{TrmtLayout}{a data frame containing the treatment layout of the experiment, if \code{genLayout} is \code{FALSE} or a list whose length is equal to the 
#'      number of trials containing the treatment layout for each trial, if \code{genLayout} is \code{TRUE}} 
#' \item{PlotNumLayout}{a matrix containing the plot number information of the experiment} 
#' \item{RepLayout}{a matrix containing the replication information of the experiment} 
#' \item{RowBlockLayout}{a matrix containing the row block information of the experiment}
#' \item{ColumnBlockLayout}{a matrix containing the column block information of the experiment}
#' 
#' @examples
#' ## Generate randomization for an experiment with 20 treatment levels in Row-Column design replicated 4 times.
#' ## Each replicate will have 4 row blocks.
#' rowcol1 <- designRowColumn(generate = list(Entry = 20), numRowBlk = 4, numRep = 4, numTrial = 1, 
#'                            genLayout = FALSE)
#' 
#' ## Generate randomization and layout for an experiment with 20 treatment levels in Row-Column design replicated 
#' ## 4 times. Each replicate will have 4 row blocks. The experiment will be arrange in a 8 x 10 field while each 
#' ## replicate will be arrange in a 4 x 5.
#' 
#' rowcol2 <- designRowColumn(generate = list(Entry = 20), numRowBlk = 4, numRep = 4, numTrial = 2, 
#'                            genLayout = TRUE, numFieldRow = 8)
#' 
# ------------------------------------------------------------------------------------

designRowColumn <- function(generate, numRowBlk = 2, numRep = 2, numTrial = 1, 
                            genLayout = FALSE, numFieldRow = 2, serpentine = FALSE, 
                            topToBottom = TRUE, display = TRUE) { #UseMethod("designRowColumn") }

#designRowColumn.default <- function(generate, numRowBlk = 2, numRep = 2, numTrial = 1, 
#                                    genLayout = FALSE, numFieldRow = 2, serpentine = FALSE, 
#                                    topToBottom = TRUE, display = TRUE) {
#  
  # --- check inputs --- #
  
  if (numRep < 2) { 
    stop("ERROR: The number of replicates should be greater than or equal to 2.")} # check number of replicates, should be >= 2
  
  if (length(generate[[1]]) == 1) { 
    tempComb <- FactorList(generate) 
  } else { 
    tempComb <- generate } # create the levels of the treatment
  
  flag <- primeNumber(length(tempComb[[1]]))
  if(flag == 1) {
    stop("ERROR: The number of treatments should not be a prime number.")}
  
  if(length(tempComb[[1]]) < 9){
    stop("ERROR: The number of treatments should be greater than or equal to 9.")}
  
  if (numRowBlk <= 1 || numRowBlk >= length(tempComb[[1]])) { 
    stop("ERROR: The number of rows per replicate should not be equal to 1 or the number of treatments.") } # check if rowPerRep is greater than 1 or
  
  if (!genLayout) {
    numFieldRow <- numRowBlk }
  
  if ((numRep * length(tempComb[[1]])) > 1500) { 
    stop("ERROR: The maximum number of experimental units that can be generated is 1500.") }   # determine the total number of experimental units limitation of DiGGer package
  
  if ((length(tempComb[[1]])*numRep)%%numFieldRow != 0) { 
    stop("ERROR: The total number of plots should be divisible by the number of field rows.") }
  
  numFieldCol <- (length(tempComb[[1]])*numRep)/numFieldRow    # determine the number of field column in the experiment
  
  if (numFieldRow %% numRowBlk != 0) { 
    stop("ERROR: The total number of plots should be divisible by the number of field rows.") }
  
  numRepRow <- numFieldRow/numRowBlk                     # determine the number of rep along the length of the field layout
  
  if (!(numRepRow %in% allFactors(numRep))) { 
    stop("ERROR: The quotient of the number of field rows and number of rows in each replicate should be a factor of the number of the replicates.") }
  
  # determine the number of columns in each replicate
  if ((length(tempComb[[1]]) %% numRowBlk) != 0) { 
    stop("ERROR: The total number of treatment levels should be divisible by the number of replicates.") }
  
  numColBlk <- length(tempComb[[1]])/numRowBlk   # determine the number of columns per replicate
  
  if (numFieldCol %% numColBlk != 0) { 
    stop("ERROR: The total number of plots should be divisible by the number of field rows.") }
  
  numRepCol <- numFieldCol/numColBlk
  
  if (numRep * length(tempComb[[1]]) != numFieldRow * numFieldCol) { 
    stop("ERROR: The total number of plots cannot be accomodated in the field experiment.") }
  
  fieldbook <- NULL
  trmtLayout <- NULL
  
  for (i in (1:numTrial)) {
    tmpResult <- randomizeRCD(numTrmt = length(tempComb[[1]]), numFieldRow, numFieldCol, 
                              numRep, numRowBlk, numColBlk,
                              trmtList = tempComb[[1]]) 
    
    # fieldbook <- rbind(fieldbook, data.frame(Trial = i, tmpResult$fieldbook))
    # trmtLayout <- rbind(trmtLayout, data.frame(Trial = i, tmpResult$plan))
    
    fieldbook <- rbind(fieldbook, data.frame(Trial = i, tmpResult$fieldbook))
    
    if (!genLayout) {
      trmtLayout <- rbind(trmtLayout, data.frame(Trial = i, tmpResult$plan))
    } else {
      if (i == 1) trmtLayout <- list()
      trmtLayout[[i]] <- tmpResult$plan
    }
    
  } ## end stmt -- for (i in (1:numTrial))
  
  fieldbook$RowBlock <- fieldbook$ROW %% numRowBlk
  fieldbook$ColumnBlock <- fieldbook$RANGE %% numColBlk
  fieldbook[fieldbook$RowBlock == 0, "RowBlock"] <- numRowBlk
  fieldbook[fieldbook$ColumnBlock == 0, "ColumnBlock"] <- numColBlk
  
  if (!genLayout) { 
    # fieldbook$PlotNumber <- as.numeric(paste(fieldbook$REP, paste(rep(0,nchar(numTrmt)), collapse=""), sep = "")) + 1:numTrmt 
    fieldbook$PlotNumber <- as.numeric(paste(fieldbook$REP, paste(rep(0,nchar(length(tempComb[[1]]))), collapse=""), sep = "")) + 1:length(tempComb[[1]]) 
    fieldbook <- fieldbook[,c("Trial", "REP", "RowBlock", "ColumnBlock","ID", "PlotNumber")]
    names(fieldbook) <- c("Trial", "Rep", "RowBlock", "ColumnBlock", names(tempComb)[1], "PlotNumber")
    tmp <- trmtLayout
    names(tmp)[2:ncol(tmp)] <- paste(rep(paste("Rep", 1:numRep, sep = ""), each = numColBlk),
                                     paste("ColumnBlock", 1:numColBlk, sep = ""), sep = "-")
    startIndex <- 2; endIndex <- numColBlk+1
    
    for (i in 1:numRepCol){
      if (i == 1){ 
        tmpVarying <- paste("c(",startIndex,":",endIndex,")", sep = "")
      } else { 
        tmpVarying <- paste(tmpVarying, ", c(", startIndex,":",endIndex,")", sep = "") }
      
      startIndex <- startIndex + numColBlk
      endIndex <- endIndex + numColBlk
    }
    
    tmp$RowBlock <- rep(1:numRowBlk, times = numTrial)
    command <- paste("reshape(tmp, idvar = c('Trial','RowBlock'), v.names = paste('Rep', 1:numRep, sep =''), timevar = 'ColumnBlock', varying = list(",tmpVarying,"),direction = 'long', new.row.names = 1:(length(tempComb[[1]])*numTrial))", sep = "")
    trmtLayout <- eval(parse(text = command))
    trmtLayout <- trmtLayout[order(trmtLayout$Trial),]
    colBlkLayout <- rowBlkLayout <- trmtLayout
    colBlkLayout[,paste("Rep", 1:numRep, sep = "")] <- trmtLayout$ColumnBlock
    rowBlkLayout[,paste("Rep", 1:numRep, sep = "")] <- trmtLayout$RowBlock
    colBlkLayout <- colBlkLayout[,c("Trial", paste("Rep", 1:numRep, sep = ""))]
    rowBlkLayout <- rowBlkLayout[,c("Trial", paste("Rep", 1:numRep, sep = ""))]
    trmtLayout <- trmtLayout[,c("Trial", paste("Rep", 1:numRep, sep = ""))]
    rownames(colBlkLayout) = rownames(rowBlkLayout) = rownames(trmtLayout) <- 1:nrow(trmtLayout)
    exptLayout <- list()
    exptLayout[[1]] <- trmtLayout
    exptLayout[[2]] <- rowBlkLayout
    exptLayout[[3]] <- colBlkLayout
    names(exptLayout) <- c("TrmtLayout", "RowBlockLayout","ColumnBlockLayout")
  } else {
    tmpResult2 <- layoutIBD(fieldbook, trmtLayout, 
                            numFieldRow, numFieldCol,
                            numRowPerRep = numRowBlk, 
                            numColPerRep = numColBlk, 
                            numBlk = NULL, 
                            numRowPerBlk = NULL, 
                            numColPerBlk = NULL,
                            serpentine = serpentine,
                            topToBottom = topToBottom)
    
    tmpResult2$fieldbook <- tmpResult2$fieldbook[,c("Trial", "REP", "RowBlock", "ColumnBlock","ID", "PlotNumber", "ROW","RANGE")]
    names(tmpResult2$fieldbook) <- c("Trial", "Rep", "RowBlock", "ColumnBlock", names(tempComb)[1], "PlotNumber", "FieldRow", "FieldColumn")
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
    
  }
  
  if (display) {
    cat(toupper("Design Properties:"),"\n",sep = "")
    cat("\t","Incomplete Block Design","\n",sep = "") 
    cat("\t","Row-Column Design","\n\n",sep = "") 
    cat(toupper("Design Parameters:"),"\n",sep = "")
    cat("\t","Number of Trials = ", numTrial, "\n",sep = "")
    cat("\t","Number of Treatments = ", length(tempComb[[1]]), "\n",sep = "")
    cat("\t","Number of Replicates = ", numRep, "\n",sep = "")
    cat("\t","Number of Rows per Replicate = ", numRowBlk, "\n",sep = "")
    cat("\t","Number of Columns per Replicate = ", numColBlk, "\n\n",sep = "")
    
    if (genLayout) {
      cat("\t","Number of Field Rows = ", numFieldRow, "\n",sep = "")
      cat("\t","Number of Field Columns = ", numFieldCol, "\n\n",sep = "")    
    }
  } ## end stmt -- if (display)
  
  if (!genLayout) { 
    return(invisible(list(fieldbook = fieldbook, plan = exptLayout)))  
  } else { 
    return(invisible(list(fieldbook = tmpResult2$fieldbook, plan = tmpResult2$plan))) }
  
} ## end stmt -- designRowColumn function

