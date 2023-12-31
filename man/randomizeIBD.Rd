% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomizeIBD.R
\name{randomizeIBD}
\alias{randomizeIBD}
\title{Randomization for Incomplete Block designs}
\usage{
randomizeIBD(
  numTrmt,
  numFieldRow,
  numFieldCol,
  numRep,
  numRowPerRep,
  numColPerRep,
  numBlk = NULL,
  blksize = NULL,
  numRowPerBlk = NULL,
  numColPerBlk = NULL,
  trmtList = NULL
)
}
\arguments{
\item{numTrmt}{number of treatment}

\item{numFieldRow}{number of field rows}

\item{numFieldCol}{number of field columns}

\item{numRep}{number of replicates}

\item{numRowPerRep}{number of row per replicate}

\item{numColPerRep}{number of column per replicate}

\item{numBlk}{number of blocks per replicate}

\item{blksize}{number of plots per block}

\item{numRowPerBlk}{number of row per block}

\item{numColPerBlk}{number of column per block}

\item{trmtList}{NULL or vector or character containing the levels of the treatment}
}
\value{
A list containing dataframe and statistical design array.
}
\description{
Generate randomization.
}
\examples{
randomizeIBD(numTrmt = 16, numRep = 3, numRowPerRep = 16, numColPerRep = 1, 
             numBlk = 4, blksize = 4, numRowPerBlk = 4, numColPerBlk = 1, 
             numFieldRow = 16, numFieldCol = 3) 

}
