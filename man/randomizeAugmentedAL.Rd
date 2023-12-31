% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomizeAugmentedAL.R
\name{randomizeAugmentedAL}
\alias{randomizeAugmentedAL}
\title{Randomization for Augmented Alpha Lattice designs}
\usage{
randomizeAugmentedAL(
  numCheck,
  numTest,
  numBlk = 2,
  numRep = 2,
  genLayout = FALSE,
  numRowPerBlk = 1,
  numColPerBlk,
  numRowPerRep = 1,
  numColPerRep,
  numFieldRow = 1,
  numFieldCol,
  serpentine = FALSE,
  topToBottom = TRUE,
  checkTrmtList = NULL,
  testTrmtList = NULL
)
}
\arguments{
\item{numCheck}{number of replicated treatments}

\item{numTest}{number of unreplicated treatments}

\item{numBlk}{number of blocks per replicate}

\item{numRep}{number of replicates}

\item{genLayout}{logical, whether a layout of the design will be generated or not}

\item{numRowPerBlk}{number of rows per block, if genLayout is TRUE}

\item{numColPerBlk}{number of columns per block, if genLayout is TRUE}

\item{numRowPerRep}{number of rows per replicate, if genLayout is TRUE}

\item{numColPerRep}{number of columns per replicate, if genLayout is TRUE}

\item{numFieldRow}{number of field rows, if genLayout is TRUE}

\item{numFieldCol}{number of field columns, if genLayout is TRUE}

\item{serpentine}{a logical variable indicating whether plot number will be arranged as serpentine order}

\item{topToBottom}{a logical variable indicating whether plot number will be written from to to bottom}

\item{checkTrmtList}{NULL or a character vector indicating the names of the replicated treatment}

\item{testTrmtList}{NULL or a character vector indicating the names of the unreplicated treatment}
}
\value{
A list containing dataframe and statistical design array.
}
\description{
Generate randomization.
}
\examples{
randomizeAugmentedAL <- function(numCheck = 6, numTest = 12, numBlk = 3, numRep = 2,  
                                 genLayout = FALSE,
                                 checkTrmtList = paste("Check",1:6,sep=""), 
                                 testTrmtList = paste("Test",1:12,sep=""))

}
