% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomizeAlphaLatticeWithDiagCheck.R
\name{randomizeAlphaLatticeWithDiagCheck}
\alias{randomizeAlphaLatticeWithDiagCheck}
\title{Randomization for Alpha Lattice design with diagonal checks}
\usage{
randomizeAlphaLatticeWithDiagCheck(
  numTrmt,
  numSpatCheck,
  numFieldRow,
  numFieldCol,
  numRep,
  numRowPerRep,
  numColPerRep,
  numBlk,
  blksize,
  numRowPerBlk,
  numColPerBlk,
  topToBottom = NULL,
  serpentine = NULL,
  trmtList = NULL
)
}
\arguments{
\item{numTrmt}{number of treatment}

\item{numSpatCheck}{number of spatial checks}

\item{numFieldRow}{number of field rows}

\item{numFieldCol}{number of field columns}

\item{numRep}{number of replicates}

\item{numRowPerRep}{number of row per replicate}

\item{numColPerRep}{number of column per replicate}

\item{numBlk}{number of blocks per replicate}

\item{blksize}{number of plots per block}

\item{numRowPerBlk}{number of row per block}

\item{numColPerBlk}{number of column per block}

\item{topToBottom}{a logical variable indicating whether plot number will be written from to to bottom}

\item{serpentine}{a logical variable indicating whether plot number will be arranged as serpentine order}

\item{trmtList}{NULL or vector or character containing the levels of the treatment}
}
\value{
A list containing dataframe and statistical design array.
}
\description{
Generate randomization.
}
\examples{
randomizeAlphaLatticeWithDiagCheck(numTrmt = 27, numSpatCheck = 1, 
                                   numFieldRow = 9, numFieldCol = 6, 
                                   numRep = 2, numRowPerRep = 9, numColPerRep = 3, 
                                   numBlk = 9, blksize = 3, numRowPerBlk = 3, numColPerBlk = 1,
                                   topToBottom = TRUE, serpentine = TRUE, 
                                   trmtList = c(paste("Entry",1:27,sep=""),paste("SpatialCheck",1,sep="")))

}
