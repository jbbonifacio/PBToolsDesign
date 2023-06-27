# -------------------------------------------------------------------------------------
# Name             : GenerateFactor
# Description      : Generate all treatment combination of several factors.
# R Version        : 4.0.1
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 2.0.0
# Modified by      : Alaine A. Gulles 
# Date Modified    : 2020.09.28
# Changes made     : 
# -------------------------------------------------------------------------------------

#' @name GenerateFactor
#' @aliases GenerateFactor
#  @aliases GenerateFactor.default
#' @title Create list of treatments or treatment combinations
#'
#' @description Generate all treatment combination of several factors.
#'
#' @param generate a list
#' @param times number of times the levels will be repeated
#'
#' @return a dataframe.
#'
#' @examples
#' GenerateFactor(generate = list(Variety = 4))
#'
# -------------------------------------------------------------------------------------

GenerateFactor <- function(generate, times = 1) { #UseMethod("GenerateFactor") }

# GenerateFactor.default <- function(generate, times = 1) {
	if (!is.list(generate)) stop("The argument 'generate' must be a list.")
	numFactor <- length(generate)
	trmtNumLevel <- NULL

	generate <- FactorList(generate)

	if (length(generate) == 1) { 
	  trmt <- data.frame(gl(length(generate[[1]]), times, labels = generate[[1]]))
	} else {
		for (i in (1:length(generate))) { 
		  trmtNumLevel <- c(trmtNumLevel, length(generate[[i]])) 
		}
		
	  trmt <- data.frame(gl(length(generate[[1]]), prod(trmtNumLevel[2:length(trmtNumLevel)])*times, labels = generate[[1]]))
		
	  for (i in (2:length(trmtNumLevel))) {
			if (i == length(trmtNumLevel)) { 
			  trmt <- cbind(trmt, data.frame(gl(trmtNumLevel[i], times, prod(trmtNumLevel[1:length(trmtNumLevel)])*times, labels = generate[[i]])))
			} else { 
			  trmt <- cbind(trmt, data.frame(gl(trmtNumLevel[i], prod(trmtNumLevel[(i+1):length(trmtNumLevel)])*times, prod(trmtNumLevel[1:length(trmtNumLevel)])*times, labels = generate[[i]]))) 
			}
	  }
	  
	}

	names(trmt) <- names(generate)
	return(trmt)
	
} # end of GenerateFactor function
