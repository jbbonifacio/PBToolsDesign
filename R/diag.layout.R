# ------------------------------------------------------------------------------------------
# Description: This function is from the ebsRtools package which returns a 
# diagonal layout, it was modified to change starting position of the diagonals
# Source: ebsRtools package version 0.2.0
#
# ------------------------------------------------------------------------------------------
#' @name diag.layout
#' @aliases diag.layout
#' @title diag.layout
#' 
#' @description Function that returns a diagonal layout
#' @param nCol Number of columns in the field, columns are the X coordinates
#' @param nRow Number of rows in the field, rows are the Y coordinates
#' @param percentChk Percentage of the total plots that should be check plots. Values must be between 0 and 0.5
#' @param plot TRUE or FALSE. If should display the field layout in a chart. Default FALSE
#' @source ebsRtools package version 0.2.0
#' 
# ------------------------------------------------------------------------------------------

diag.layout <- function (nCol, nRow, percentChk, plot = FALSE) { #UseMethod("diag.Layout") }
  
  #diag.layout.default <- function (nCol, nRow, percentChk, plot = FALSE) {
  k <- as.integer(1/percentChk)
  if (percentChk %in% c(0.08, 0.11, 0.14, 0.15)) {
    if (percentChk == 0.11) {
      l <- 8
    }
    if (percentChk == 0.14) {
      l <- 3
    }
    if (percentChk %in% c(0.08, 0.15)) {
      l <- 6
    }
  } else {
    l <- 4
  }
  
  startPos <- sample(k,1)
  pattern <- rep(0, k)
  pattern[startPos] <- 1
  
  tmpnCol <- length(rep(pattern, (nCol/length(pattern)) + 1))
  layout <- matrix(NA, nRow, (tmpnCol))
  layout[1, ] <- rep(pattern, (nCol/length(pattern)) + 1)
  
  for (i in c(2:nRow)) {
    layout[i, ] <- layout[i - 1, ][c(c(l:tmpnCol, c(1:(l - 1))))]
    if (i == nRow) {
      layout <- layout[, 1:nCol]
    }
  }
  
  layout <- as.data.frame(layout)
  layout <- layout[with(layout, order(as.numeric(rownames(layout)), 
                                      decreasing = F)), ]
  layout <- as.matrix(layout)
  colnames(layout) <- c(1:nCol)
  realPrct <- sum(layout == 1)/(sum(layout == 1) + sum(layout == 0))
  if (plot) {
    heatmap(layout, Rowv = NA, Colv = NA, scale = "none", 
            xlab = "field col", ylab = "field row", 
            main = paste("check plots", realPrct * 100, 
                         "%"), col = c("white", "black"), 
            add.expr = abline(h = 0.5 + c(1:nRow), v = 0.5 + 
                                c(1:nCol)), margins = c(3, 3))
  }
  return(layout)
}
