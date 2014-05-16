\name{get.cu.weights}
\alias{get.cu.weights}
\title{Extracts weights for Cu pups
: from various tables in the cutagnew.mdb ACCESS
database and returns a dataframe with ancillary fields added.}
\usage{
get.cu.weights(fdir = NULL, ENYears = c(1976, 1983, 1984, 1987, 1992, 1997,
  1998))
}
\arguments{
  \item{fdir}{directory for cutagnew.mdb}

  \item{ENYears}{values of years that will be flagged as
  ENSO years}
}
\value{
dataframe of weights for Cu
}
\description{
Extracts weights for Cu pups : from various tables in the
cutagnew.mdb ACCESS database and returns a dataframe with
ancillary fields added.
}
\examples{
fdir=system.file(package="CIPinnipedAnalysis")
source(file.path(fdir,"CuWeightAnalysis.r"))
}
\author{
Jeff Laake
}

