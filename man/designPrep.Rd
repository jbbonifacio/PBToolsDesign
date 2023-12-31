% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/designPRepDesign.R
\name{designPrep}
\alias{designPrep}
\title{P-Rep Design}
\usage{
designPrep(
  numTrmtPerTrmtGrp,
  numRepPerTrmtGrp,
  numTrmtPerGrp,
  numTrial = 1,
  genLayout = FALSE,
  numFieldRow = 1,
  treatName = NULL,
  listName = NULL,
  serpentine = FALSE,
  topToBottom = TRUE,
  display = TRUE
)
}
\arguments{
\item{numTrmtPerTrmtGrp}{number of treatments per replicate group}

\item{numRepPerTrmtGrp}{number of replicates per treatment group}

\item{numTrmtPerGrp}{number of treatments per group}

\item{numTrial}{number of trials (randomization set-ups with the complete set of entries)}

\item{genLayout}{logical, whether a layout of the design will be generated or not}

\item{numFieldRow}{number of field rows}

\item{treatName}{NULL or vector of names associated with the treatment}

\item{listName}{column name associated with the list provided (treatName)}

\item{serpentine}{logical, whether plot number will be arranged as serpentine order}

\item{topToBottom}{logical, whether plot number will be written from to to bottom}

\item{display}{logical, whether randomization parameters will be displayed}
}
\value{
A list containing the following components:
\item{fieldbook}{a data frame}
\item{plan}{a data frame, if \code{genLayout} is \code{FALSE} or a list containing the following components 
if \code{genLayout} is \code{TRUE}}
\item{TrmtLayout}{a list whose length is equal to the 
     number of trials containing the treatment layout for each trial} 
\item{PlotNumLayout}{a matrix containing the plot number layout of the experiment}
}
\description{
Generate randomization and layout for P-Rep Design.
}
\details{
If \code{genLayout} is \code{TRUE}, then parameters \code{numFieldRow} and \code{numFieldCol} must be specified.  
Values of \code{numFieldRow} and \code{numFieldCol} should be a factor of the total number of plots.
}
\examples{
## Generate randomization for an experiment with 150 treatment levels in P-rep design.
## There will be 2 groups (test and check)
## Test group 1 - 175 entries replicated 1 time
## Check group 1 - 5 entries replicated 5 times
prep1 <- designPrep(numRepPerTrmtGrp = c(1,5), numTrmtPerTrmtGrp = c(175,5), numTrmtPerGrp = c(175,5), numTrial = 1,
                   genLayout = FALSE, treatName = NULL, listName = "Entry",
                   display = TRUE)
                   
## The experiment will be arranged in a 20 x 10 field in Left to Right Serpentine order.
prep2 <- designPrep(numRepPerTrmtGrp = c(1,2,5), numTrmtPerTrmtGrp = c(115,30,5), numTrmtPerGrp = c(145,5), numTrial = 1,
                   genLayout = TRUE, numFieldRow = 20, treatName = NULL, listName = "Entry",
                   serpentine = TRUE, topToBottom = FALSE, display = TRUE)
}
