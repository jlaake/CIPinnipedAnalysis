\name{do.early.pup.mortality}
\alias{do.early.pup.mortality}
\title{Early pup mortality estimation
For Zc and Cu, creates data tables in the ACCESS database with the early pup
mortality estimates during the season in each year.  Also constructs pdf
plots of those values.}
\usage{
  do.early.pup.mortality(fdir = "")
}
\arguments{
  \item{fdir}{directory for CIPinnipedCensusQuery.mdb}
}
\value{
  None
}
\description{
  Uses the \code{\link{mortality.stats}}function to create
  two data tables of early pup mortality for Cu and Zc in
  CIPinnipedCensusQuery.mdb and it creates plots in pdf
  files. If there are any errors then the mortality tables
  will not be created. If successful the pdf files are
  created in this directory which can then be copied to
  Presentations/PinnipedPlots.
}
\author{
  Jeff Laake
}

