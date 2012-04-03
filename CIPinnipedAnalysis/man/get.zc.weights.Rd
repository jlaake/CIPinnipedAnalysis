\name{get.zc.weights}
\alias{get.zc.weights}
\title{Extracts weights for Zc pups
:from various tables in the BrandMaster.mdb ACCESS database and returns a dataframe with ancillary fields added.}
\usage{
  get.zc.weights(fdir = NULL,
    ENYears = c(1976, 1983, 1984, 1987, 1992, 1997, 1998))
}
\arguments{
  \item{fdir}{directory for BrandMaster.mdb}

  \item{ENYears}{values of years that will be flagged as
  ENSO years}
}
\value{
  dataframe of weights for Zc
}
\description{
  Extracts weights for Zc pups :from various tables in the
  BrandMaster.mdb ACCESS database and returns a dataframe
  with ancillary fields added.
}
\examples{
fdir=system.file(package="CIPinnipedAnalysis")
source(file.path(fdir,"ZcWeightAnalysis.r"))
}
\author{
  Jeff Laake
}

