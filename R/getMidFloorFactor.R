# -------------------------------------------------------------------------------------
# Name             : getMidFloorFactor 
# Description      : Get the lower middle factor of a number
# R Version        : 4.1.0 
# -------------------------------------------------------------------------------------
# Author           : Justine B. Dayrit 
# Author Email     : j.bonifacio@irri.org
# Script Version   : 1.0.0
# Modified by      : 
# Date Modified    : 
# Changes made     : 
# -------------------------------------------------------------------------------------
#' @name getMidFloorFactor
#' @aliases getMidFloorFactor
#  @aliases getMidFloorFactor.default
#' @title getMidFloorFactor
#'
#' @description Get the lower middle factor of a number
#'
#' @param x number 
#' 
#' @return an integer
#' 
#' @examples
#' getMidFloorFactor(20)
# -------------------------------------------------------------------------------------

getMidFloorFactor <- function(x) { #UseMethod("getMidFloorFactor") }

#getMidFloorFactor.default <- function(x) {
  factors <- NULL
  for(i in 1:x) {
    if((x %% i) == 0) {
      factors <- c(factors,i)
    }
  }
  
  if(length(factors) %% 2 == 0){
    return(as.numeric(factors[length(factors)/2]))
  } else{
    return(as.numeric(median(factors)))
  }
  
} # end of getMidFloorFactor function