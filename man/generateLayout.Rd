% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateLayout.R
\name{generateLayout}
\alias{generateLayout}
\title{Layout Creation}
\usage{
generateLayout(
  fieldbook,
  numFieldRow = 1,
  numRowPerBlk = 1,
  serpentine = FALSE,
  topToBottom = TRUE
)
}
\arguments{
\item{fieldbook}{a dataframe which is a result from the function \code{randomizeRCBD}}

\item{numFieldRow}{number of field rows}

\item{numRowPerBlk}{number of rows per block (or replicate)}

\item{serpentine}{a logical variable indicating whether plot number will be arranged as serpentine order}

\item{topToBottom}{a logical variable indicating whether plot number will be written from to to bottom}
}
\value{
A list containing the following components:
\item{fieldbook}{a data frame}
\item{plan}{a list containing the following components:}
\item{TrmtLayout}{a list whose length is equal to the number of trials containing the treatment layout for each trial}
\item{PlotNumLayout}{a matrix containing the plot number layout of the experiment}
\item{BlockLayout}{a matrix containing the replication layout of the experiment}
}
\description{
Creates layout from an experimental design in Randomized Complete Block.
}
\examples{
myFieldbook <- randomizeRCBD(generate = list(Entry = 1:10), numRep = 4, numTrial = 1)
myLayout <- generateLayout(myFieldbook, numFieldRow = 10, numRowPerBlk = 5)

}
