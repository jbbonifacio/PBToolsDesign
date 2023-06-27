# -------------------------------------------------------------------------------------
# Name             : randomizeRCBD
# Description      : Generate randomization for randomized complete block design (RCBD)
#                    for single factor or factorial experiments.
# R Version        : 4.0.1 
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 2.0.1
# Modified by      : Justine B. Bonifacio
# Date Modified    : 2020.11.11
# Changes made     : converted Block from factor to numeric
# -------------------------------------------------------------------------------------
#' @name randomizeRCBD
#' @aliases randomizeRCBD
#  @aliases randomizeRCBD.default
#' @title Randomization for Randomized Complete Block Design (RCBD)
#'
#' @description Generate randomization for randomized complete block design (RCBD)
#'             for single factor or factorial experiments.
#'
#' @param generate list of entries to be randomized
#' @param numBlk number of replicates or blocks
#' @param numTrial number of trials (randomization set-ups with the complete set of entries)
#'
#' @return A dataframe.
#'
#' @examples
#' randomizeRCBD(generate = list(Variety = 4), numBlk = 2, numTrial = 1)
#'
# -------------------------------------------------------------------------------------

randomizeRCBD <- function(generate, numBlk = 2, numTrial = 1) { #UseMethod("randomizeRCBD") }

# randomizeRCBD.default <- function(generate, numBlk = 2, numTrial = 1) {

  if (is.null(numTrial) || numTrial < 1 || is.character(numTrial) || length(numTrial) > 1) { stop("The argument 'numTrial' should be a single value greater than or equal to 1.") }
  if (is.null(numBlk) || numBlk < 2 || is.character(numBlk) || length(numBlk) > 1) { stop("The argument 'numBlk' should be a single value greater than or equal to 2.") }
  if (missing(generate)) { stop("The argument 'generate' is missing.") }
  if (!is.list(generate)) { stop("The argument 'generate' must be a list.") }

  tempComb <- GenerateFactor(generate, times = 1)
  randomize <- NULL

  for (i in (1:numTrial)) {
    for (j in (1:numBlk)) {
      temp <- data.frame(Trial = as.character(i), Block = as.character(j), tempComb, tempPlotNum = sample(nrow(tempComb), nrow(tempComb), replace = FALSE))
      temp <- temp[order(temp[,"tempPlotNum"]),]
      plotLabel <- as.numeric(paste(j, paste(c(rep(0, max(nchar(1:nrow(tempComb))))), collapse = ""), sep = ""))+1:nrow(tempComb)
      if (ncol(tempComb) > 1) { trmtLabel <- eval(parse(text = paste("paste(temp[,'", paste(names(temp)[3:(ncol(temp)-1)], collapse = "'],' ',temp[,'", sep = ""),"'], sep = '')", sep = "")))
      } else { trmtLabel <- temp[,3] }
      randomize <- rbind(randomize, cbind(temp, trmtLabel, plotLabel))
    }
  }
  
  # randomize <- randomize[,-I(c(match(c("trmtLabel", "tempPlotNum"), names(randomize))))]
  randomize <- randomize[,-I(c(match(c("tempPlotNum"), names(randomize))))]
  names(randomize)[ncol(randomize)] <- c("PlotNumber")
  randomize$Block <- as.numeric(randomize$Block)
  randomize <- randomize[order(randomize$Trial, randomize$Block, randomize$PlotNumber),]
  rownames(randomize) <- 1:nrow(randomize)

  return(fieldbook = randomize)
}
