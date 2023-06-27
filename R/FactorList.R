# -------------------------------------------------------------------------------------
# Name             : FactorList
# Description      : Generate all treatment combination of several 
# R Version        : 4.0.1
# -------------------------------------------------------------------------------------
# Author           : Alaine A. Gulles 
# Author Email     : a.gulles@irri.org
# Script Version   : 2.0.0
# Modified by      : Alaine A. Gulles
# Date Modified    : 2020.09.28
# Changes made     : 
# -------------------------------------------------------------------------------------
#' @name FactorList
#' @aliases FactorList
#  @aliases FactorList.default
#' @title Create levels of a factor
#'
#' @description Generate all treatment combination of several factors.
#'
#' @param generate a list containing the different levels of treatment combination
#'
#' @return a dataframe.
#'
#' @examples
#' FactorList(generate = list(Entry = 1:10))
#'
# -------------------------------------------------------------------------------------

FactorList <- function(generate) { #UseMethod("FactorList") }

# FactorList.default <- function(generate){
	if (!is.list(generate)) stop("The argument 'generate' must be a list.")
  
	for (i in (1:length(generate))) {
		if (length(generate[[i]]) == 1 && is.numeric(generate[[i]])) {
			generate[[i]] <- paste(names(generate)[i],1:generate[[i]], sep = "")
		}
	}
  
	return(generate)
} # end of FactorList function
