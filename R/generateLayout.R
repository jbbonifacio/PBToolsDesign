# -------------------------------------------------------------------------------------
# Name             : generateLayout
# Description      : Generate layout
# R Version        : 4.0.1
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 3.0.0
# Modified by      : Justine B. Bonifacio
# Date Modified    : 2020.09.28
# Changes made     : 
# -------------------------------------------------------------------------------------
#' @name generateLayout
#' @aliases generateLayout
#  @aliases generateLayout.default
#' @title Layout Creation
#'
#' @description Creates layout from an experimental design in Randomized Complete Block.
#'
#' @param fieldbook a dataframe which is a result from the function \code{randomizeRCBD}
#' @param numFieldRow number of field rows
#' @param numRowPerBlk number of rows per block (or replicate)
#' @param serpentine a logical variable indicating whether plot number will be arranged as serpentine order
#' @param topToBottom a logical variable indicating whether plot number will be written from to to bottom
#'
#' @return A list containing the following components:
#' \item{fieldbook}{a data frame}
#' \item{plan}{a list containing the following components:}
#' \item{TrmtLayout}{a list whose length is equal to the number of trials containing the treatment layout for each trial}
#' \item{PlotNumLayout}{a matrix containing the plot number layout of the experiment}
#' \item{BlockLayout}{a matrix containing the replication layout of the experiment}
#'
#' @examples
#' myFieldbook <- randomizeRCBD(generate = list(Entry = 1:10), numRep = 4, numTrial = 1)
#' myLayout <- generateLayout(myFieldbook, numFieldRow = 10, numRowPerBlk = 5)
#'
# -------------------------------------------------------------------------------------

generateLayout <- function(fieldbook, numFieldRow = 1, numRowPerBlk = 1,
                           serpentine = FALSE, topToBottom = TRUE) { #UseMethod("generateLayout") }

# generateLayout.default <- function(fieldbook, numFieldRow = 1, numRowPerBlk = 1,
#                                    serpentine = FALSE, topToBottom = TRUE) {
  
  fieldbook[,1] <- as.factor(fieldbook[,1])
  fieldbook[,2] <- as.factor(fieldbook[,2])
  fieldbook[,"trmtLabel"] <- as.factor(fieldbook[,"trmtLabel"])

  numRepRow <- numFieldRow/numRowPerBlk
  numRepCol <- nlevels(fieldbook[,2])/numRepRow
  colPerRep <- (nlevels(fieldbook[,"trmtLabel"])/numRowPerBlk)
  numFieldCol <- (nlevels(fieldbook[,"trmtLabel"])*nlevels(fieldbook[,2]))/numFieldRow
  
  trmtLayout <- list()
  newFieldbook <- NULL
  plotnumLayout <- matrix(0, nrow = numFieldRow,
                          ncol = (nlevels(fieldbook[,"trmtLabel"])*nlevels(fieldbook[,2]))/numFieldRow)
  repLayout <- matrix(0, nrow = numFieldRow,
                      ncol = (nlevels(fieldbook[,"trmtLabel"])*nlevels(fieldbook[,2]))/numFieldRow)
  
  rowIndexLL <- 0
  rowIndexUL <- 0
  colIndexLL <- 0
  colIndexUL <- 0
  
  for (i in 1:nlevels(fieldbook[,1])) {
    temp <- fieldbook[fieldbook[1] == levels(fieldbook[,1])[i],]
    trmtLayout[[i]] <- matrix(0, nrow = numFieldRow,
                              ncol = (nlevels(fieldbook[,"trmtLabel"])*nlevels(fieldbook[,2]))/numFieldRow)

    for (j in 1:nlevels(temp[,2])) {
      tempRepInfo <- temp[temp[,2] == levels(temp[,2])[j],]
      if (i == 1) tempRep <- matrix(rep(j, length(tempRepInfo[,"trmtLabel"])), nrow = numRowPerBlk,
                                    ncol = colPerRep)

      # orientation of the plotnumber topToBottom or leftToRight
      if (topToBottom) {
        tempOrder <- serpentine
        tempTrmt <- matrix(tempRepInfo[,"trmtLabel"], nrow = numRowPerBlk, ncol = colPerRep, byrow = FALSE)
        
        if (i == 1) tempPlotnum <- matrix(tempRepInfo[,"PlotNumber"], nrow = numRowPerBlk, ncol = colPerRep, byrow = FALSE)
        
        if (serpentine & colPerRep != 1) { 
          for (k in seq(2, colPerRep, by = 2)) {
            tempTrmt[, k] <- rev(tempTrmt[, k])
            if (i == 1) tempPlotnum[, k] <- rev(tempPlotnum[, k])
          }
        }
        if(serpentine & colPerRep == 1){
          tempOrder <- FALSE
        }
      } else {
        tempOrder <- serpentine
        tempTrmt <- matrix(tempRepInfo[,"trmtLabel"], nrow = numRowPerBlk, ncol = colPerRep, byrow = TRUE)
        
        if (i == 1) tempPlotnum <- matrix(tempRepInfo[,"PlotNumber"], nrow = numRowPerBlk, ncol = colPerRep, byrow = TRUE)
        
        if (serpentine & numRowPerBlk != 1) {
          for (k in seq(2, numRowPerBlk, by = 2)) {
            tempTrmt[k, ] <- rev(tempTrmt[k, ])
            if (i == 1) tempPlotnum[k, ] <- rev(tempPlotnum[k, ])
          }
        }
        if(serpentine & numRowPerBlk == 1){
          tempOrder <- FALSE
        }
      }
      
      range <- getRange(index = j, topToBottom = topToBottom, serpentine = tempOrder,
                        numRowPerGrp = numRowPerBlk, numColPerGrp = colPerRep,
                        numGrpRow = numRepRow, numGrpCol = numRepCol,
                        numRow = numFieldRow, numCol = numFieldCol,
                        rowIndexLL = rowIndexLL, rowIndexUL = rowIndexUL,
                        colIndexLL = colIndexLL, colIndexUL = colIndexUL)
      
      rowIndexLL <- range$rowIndexLL
      rowIndexUL <- range$rowIndexUL
      colIndexLL <- range$colIndexLL
      colIndexUL <- range$colIndexUL
      
      trmtLayout[[i]][rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tempTrmt
      
      if (i == 1) {
        plotnumLayout[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tempPlotnum
        repLayout[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <-  tempRep
      }

    } # end of for j stmt

    if(topToBottom){
      if(serpentine & colPerRep == 1 & numFieldCol != 1){
        for(m in seq(2, numFieldCol, by = 2)){
          trmtLayout[[i]][ ,m] <- rev(trmtLayout[[i]][ ,m])
          if (i == 1) {
            plotnumLayout[ ,m] <- rev(plotnumLayout[ ,m])
            repLayout[ ,m] <- rev(repLayout[ ,m])
          }
        }
      }
    } else{
      if(serpentine & numRowPerBlk == 1 & numFieldRow != 1){
        for(m in seq(2, numFieldRow, by = 2)){
          trmtLayout[[i]][m, ] <- rev(trmtLayout[[i]][m, ])
          if (i == 1) {
            plotnumLayout[m, ] <- rev(plotnumLayout[m, ])
            repLayout[m, ] <- rev(repLayout[m, ])
          }
        }
      }
    }
    
    tempFB <- merge(as.data.frame.table(repLayout), as.data.frame.table(trmtLayout[[i]]), by = c("Var1", "Var2"))
    names(tempFB)[3:4] <- c(names(temp)[2], "trmtLabel")
    tempFB <- merge(tempFB, as.data.frame.table(plotnumLayout), by = c("Var1", "Var2"))
    names(tempFB)[ncol(tempFB)] <- names(temp)[ncol(temp)]

    newFieldbook <- rbind(newFieldbook, merge(temp, tempFB))
    dimnames(trmtLayout[[i]]) <- list(paste("FieldRow", 1:nrow(trmtLayout[[i]]),sep = ""),
                                      paste("FieldColumn", 1:ncol(trmtLayout[[i]]), sep = ""))
  } # end of for i stmt

  newFieldbook[,"Var1"] <- as.numeric(newFieldbook[,"Var1"])
  newFieldbook[,"Var2"] <- as.numeric(newFieldbook[,"Var2"])
  newFieldbook <- newFieldbook[,c(names(temp), "Var1", "Var2")]
  names(newFieldbook)[(ncol(newFieldbook)-1):ncol(newFieldbook)] <- c("FieldRow", "FieldColumn")
  newFieldbook <- newFieldbook[,-I(c(match(c("trmtLabel"), names(newFieldbook))))]
  dimnames(plotnumLayout) <- list(paste("FieldRow", 1:nrow(plotnumLayout),sep = ""),
                                  paste("FieldColumn", 1:ncol(plotnumLayout), sep = ""))
  dimnames(repLayout) <- list(paste("FieldRow", 1:nrow(repLayout),sep = ""),
                              paste("FieldColumn", 1:ncol(repLayout), sep = ""))

  exptLayout <- list()
  exptLayout[[1]] <- trmtLayout
  names(exptLayout[[1]]) <- paste("Trial", levels(fieldbook[,"Trial"]), sep = "")
  exptLayout[[2]] <- plotnumLayout
  exptLayout[[3]] <- repLayout
  names(exptLayout) <- c("TrmtLayout", "PlotNumLayout", "BlockLayout")

  return(invisible(list(fieldbook = newFieldbook, plan = exptLayout)))

} # end of generateLayout function
