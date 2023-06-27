# ------------------------------------------------------------------------------------------
# Description: This function returns 1 if the number is prime else it will return 0.
# Script from the source was slightly modified to conform with the requirement
# Source: https://www.datamentor.io/r-programming/examples/prime-number/
#
# Researched by: Justine Bonifacio
# ------------------------------------------------------------------------------------------
#' @name primeNumber
#' @aliases primeNumber
#' @title primeNumber
#' 
#' @description Function that returns 1 if the number is prime else it will return 0.
#' @param num number
#' @source https://www.datamentor.io/r-programming/examples/prime-number/
#' 
# ------------------------------------------------------------------------------------------

primeNumber <- function(num){
  
  flag <- 0
  # prime numbers are greater than 1
  if(num > 1) {
    # check for factors
    flag <- 1
    for(i in 2:(num-1)) {
      if ((num %% i) == 0) {
        flag <- 0
        break
      }
    }
  } 
  if(num == 2)    flag <- 1
  
  return(flag)
}
  