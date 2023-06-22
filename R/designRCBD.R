# -------------------------------------------------------------------------------------
# Name             : designRCBD
# Description      : Generate randomization and layout for randomized complete block
#                    design (RCBD) for single factor or factorial experiments.
# R Version        : 4.0.1 
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 4.0.1
# Modified by      : Justine B. Bonifacio
# Date Modified    : 2021.11.24
# Changes made     : added invisible function to avoid displaying of results
#                    when genLayout is TRUE
# -------------------------------------------------------------------------------------
#' @title Randomized Complete Block Design (RCBD)
#' @aliases designRCBD
#  @aliases designRCBD.default
#' @description Generate randomization and layout for randomized complete block design (RCBD)
#' for single factor or factorial experiments.
#'
#' @param generate list of entries to be randomized
#' @param numBlk number of replicates or blocks
#' @param numTrial number of trials (randomization set-ups with the complete set of entries)
#' @param genLayout logical, whether a layout of the design will be generated
#' @param numFieldRow number of field rows
#' @param numRowPerBlk number of rows per replicate or block
#' @param serpentine a logical variable indicating whether plot number will be arranged as serpentine order
#' @param topToBottom a logical variable indicating whether plot number will be written from to to bottom
#' @param display a logical variable indicating whether randomization parameters will be displayed
#'
#' @details
#' If \code{genLayout} is \code{TRUE}, then parameters \code{numFieldRow} and \code{numRowPerBlk} must be specified.
#' Values of \code{numRowPerBlk} should be a factor of \code{numFieldRow}.
#'
#' @return A list containing the following components:
#' \item{fieldbook}{a data frame}
#' \item{plan}{a data frame, if \code{genLayout} is \code{FALSE} or a list containing the following components
#' if \code{genLayout} is \code{TRUE}:}
#' \item{TrmtLayout}{a list whose length is equal to the number of trials containing the treatment layout for each trial}
#' \item{PlotNumLayout}{a matrix containing the plot number layout of the experiment}
#' \item{BlockLayout}{a matrix containing the replication layout of the experiment}
#'
#' @examples
#' ## Generate randomization for an experiment with 10 levels in RCBD replicated 4 times
#' rcbd1a <- designRCBD(generate = list(Entry = 1:10), numBlk = 4, numTrial = 1, genLayout = FALSE)
#'
#' ## Generate randomization and layout in RCBD for an experiment with l0 levels replicated 4 times in 2 trials
#' rcbd1b <- designRCBD(generate = list(Entry = 10), numBlk = 4, numTrial = 2, genLayout = TRUE,
#'    numFieldRow = 10, numRowPerBlk = 5)
#'
#' ## Generate randomization and layout in RCBD for a 5 x 4 factorial experiment replicated 4 times in 2 trials
#' varietyLevel <- paste("V", 1:5, sep = "")
#' fertLevel <- paste("F", 1:4, sep = "")
#' rcbd2a <- designRCBD(generate = list(Variety = varietyLevel, Fertilizer = fertLevel), numBlk = 4,
#'     numTrial = 2, genLayout = TRUE, numFieldRow = 10, numRowPerBlk = 5, serpentine = TRUE)
# -------------------------------------------------------------------------------------

#library(crayon)
designRCBD <- function(generate, numBlk = 2, numTrial = 1,
                       genLayout = FALSE, numFieldRow = 1, numRowPerBlk = 1,
                       serpentine = FALSE, topToBottom = TRUE,
                       display = TRUE) { #UseMethod("designRCBD") }

# designRCBD.default <- function(generate, numBlk = 2, numTrial = 1,
#                                genLayout = FALSE, numFieldRow = 1, numRowPerBlk = 1,
#                                serpentine = FALSE, topToBottom = TRUE,
#                                display = TRUE) {

  if (is.null(numTrial) || numTrial < 1 || is.character(numTrial) || length(numTrial) > 1) {
    stop("ERROR: The argument 'numTrial' should be a single value greater than or equal to 1.") }
  
  if (is.null(numBlk) || numBlk < 2 || is.character(numBlk) || length(numBlk) > 1) {
    stop("ERROR: The argument 'numBlk' should be a single value greater than or equal to 2.") }
  
  if (missing(generate)) { 
    stop("ERROR: The argument 'generate' is missing.") }
  
  if (!is.list(generate)) { 
    stop("ERROR: The argument 'generate' must be a list.") }

  tempComb <- GenerateFactor(generate, times = 1)

  if (nrow(tempComb) < 2) { 
    stop("ERROR: The number of treatments should be at least 2.") }
  
  if (genLayout) {
    if (numRowPerBlk > numFieldRow) { 
      stop("ERROR: The number of field rows should be equal to or greater than the number of rows per block.") }
    
    if (nrow(tempComb)%%numRowPerBlk != 0) { 
      stop("ERROR: The total number of plots per block should be divisible by the number of rows per block.") }
    
    if ((nrow(tempComb)*numBlk)%%numFieldRow != 0) { 
      stop("ERROR: The total number of plots should be divisible by the number of field rows.") }
    
    if (!((numFieldRow/numRowPerBlk) %in% allFactors(numBlk))) { 
      stop("ERROR: The quotient of the number of field rows and number of rows in each block should be a factor of the number of the blocks.") }
  }

  fieldbook <- randomizeRCBD(generate, numBlk, numTrial)

  if (display) {
    cat(toupper("Design Properties:"),"\n",sep = "")
    
    if (ncol(tempComb) == 1) { 
      cat("\t","Single Factor","\n",sep = "") 
    } else { 
      cat("\t","Factorial Design","\n",sep = "") }
    
    cat("\t","Randomized Complete Block Design","\n\n",sep = "")
    cat(toupper("Design Parameters:"),"\n",sep = "")
    cat("\t","Number of Trials = ", numTrial, "\n",sep = "")
    cat("\t","Number of Blocks (or Replicates) = ", numBlk, "\n",sep = "")
    
    if (ncol(tempComb) == 1) {
      cat("\t","Treatment Name = ", names(tempComb)[1], "\n",sep = "")
      cat("\t","Treatment Levels = ", sep = "")
      
      if (nlevels(tempComb[,1]) <= 5) { 
        cat(paste(levels(tempComb[,1]), collapse = ", ", sep = ""), sep = "")
      } else {
        cat(paste(levels(tempComb[,1])[1:3], collapse = ", ", sep = ""), sep = "")
        cat(paste(", ...,", levels(tempComb[,1])[nlevels(tempComb[,1])]), sep = "")
      }
      
      cat("\n\n")
    } else {
      for (i in (1:ncol(tempComb))) {
        cat("\t","Factor ",i," = ", names(tempComb)[i], "\n",sep = "")
        cat("\t","Levels = ", sep = "")
        
        if (nlevels(tempComb[,i]) <= 5) { 
          cat(paste(levels(tempComb[,i]), collapse = ", ", sep = ""), sep = "")
        } else {
          cat(paste(levels(tempComb[,i])[1:3], collapse = ", ", sep = ""), sep = "")
          cat(paste(", ...,", levels(tempComb[,i])[nlevels(tempComb[,i])]), sep = "")
        }
        
        cat("\n")
      }
      cat("\n")
    }
  }
  
  if (!genLayout) {
    fieldbook <- data.frame(index = rep(1:nlevels(fieldbook[,"trmtLabel"]), numBlk*numTrial), fieldbook)
    plan <- reshape(fieldbook, v.names = "trmtLabel", idvar = c("index","Trial"),
                    timevar = "Block", direction = "wide",
                    drop = names(fieldbook)[c(4:(match("trmtLabel", names(fieldbook))-1), ncol(fieldbook))])
    plan <- plan[2:ncol(plan)]
    names(plan)[2:ncol(plan)] <- paste("Block",1:numBlk, sep = "")
    fieldbook <- fieldbook[,-I(c(match(c("trmtLabel", "index"), names(fieldbook))))]
    
    if (nrow(tempComb)>30) { cat(green("WARNING: The total number of treatments is too large for a RCBD. You might want to consider selecting a different design.","\n\n"))}
    
    return(invisible(list(fieldbook = fieldbook, plan = plan)))
  } else {
    result <- generateLayout(fieldbook, numFieldRow, numRowPerBlk, serpentine, topToBottom)
    result$fieldbook <- result$fieldbook[order(result$fieldbook$Trial, result$fieldbook$PlotNumber),]
    rownames(result$fieldbook) <- 1:nrow(result$fieldbook)
    
    if (display) {
      cat("\t","Number of Field Row = ", numFieldRow, "\n",sep = "")
      cat("\t","Number of Field Column = ", ncol(result[[2]][[2]]), "\n\n",sep = "")
    }
    
    if (nrow(tempComb)>30) { cat(green("WARNING: The total number of treatments is too large for a RCBD. You might want to consider selecting a different design.","\n\n"))}
    
    return(invisible(result))
  }
  
} ## end stmt -- designRCBD function
