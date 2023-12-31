% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/designAlphaLatticeDiGGer.R
\name{designAlphaLatticeDiGGer}
\alias{designAlphaLatticeDiGGer}
\title{Alpha Lattice Design using DiGGer package}
\usage{
designAlphaLatticeDiGGer(
  generate,
  numBlk = 2,
  numRep = 2,
  numTrial = 1,
  genLayout = FALSE,
  numRowPerBlk = 1,
  numRowPerRep = 1,
  numFieldRow = 1,
  serpentine = FALSE,
  topToBottom = TRUE,
  display = TRUE
)
}
\arguments{
\item{generate}{list of entries to be randomized}

\item{numBlk}{number of blocks}

\item{numRep}{number of replicates}

\item{numTrial}{number of trials}

\item{genLayout}{logical, whether a layout of the design will be generated}

\item{numRowPerBlk}{number of rows per block}

\item{numRowPerRep}{number of rows per replicate}

\item{numFieldRow}{number of field rows}

\item{serpentine}{a logical variable indicating whether plot number will be arranged as serpentine order}

\item{topToBottom}{a logical variable indicating whether plot number will be written from to to bottom}

\item{display}{a logical variable indicating whether randomization parameters will be displayed}
}
\value{
A list containing the following components:
\item{fieldbook}{a data frame}
\item{plan}{a list containing the following components: TrmtLayout and BlockLayout
if \code{genLayout} is \code{FALSE} or a list containing the following components
if \code{genLayout} is \code{TRUE}:}
\item{TrmtLayout}{a data frame containing the treatment layout of the experiment, 
     if \code{genLayout} is \code{FALSE} or a list whose length is equal to the 
     number of trials containing the treatment layout for each trial, 
     if \code{genLayout} is \code{TRUE}} 
\item{PlotNumLayout}{a matrix containing the plot number layout of the experiment} 
\item{RepLayout}{a matrix containing the replication layout of the experiment} 
\item{BlockLayout}{a data frame containing the block layout of the experiment,
     if \code{genLayout} is \code{FALSE} or a matrix containing the block layout of the experiment,
     if \code{genLayout} is \code{TRUE}}
}
\description{
Generate randomization and layout for Alpha Lattice Design using 
             DiGGer package.
}
\details{
If \code{genLayout} is \code{TRUE}, then parameters \code{numFieldRow}, \code{numRowPerRep}  and \code{numRowPerBlk} 
must be specified.  Values of \code{numRowPerBlk} should be a factor of \code{numRowPerRep} while values of 
\code{numRowPerRep} should be a factor of \code{numFieldRow}.
}
\examples{
## Generate randomization for an experiment with 24 treatment levels in Alpha Lattice replicated 4 times.
## Each replicate will have 4 blocks.
alpha1 <- designAlphaLatticeDiGGer(generate = list(Entry = 24), numBlk = 4, numRep = 4, numTrial = 1, 
                                   genLayout = FALSE)

## Generate randomization and layout for an experiment with 24 levels in Alpha Lattice replicated 4 times.
## Each replicate will have 4 blocks. The experiment will be arrange in a 8 x 12 field.
## Each replicate will be arrange in a 4 x 6, while each block will be arrange in 2 x 3.
alpha2 <- designAlphaLatticeDiGGer(generate = list(Entry = 24), numBlk = 4, numRep = 4, numTrial = 1, 
                                   genLayout = TRUE, numRowPerBlk = 2, numRowPerRep = 4, numFieldRow = 8)

}
