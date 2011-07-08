\name{do.pup.production}
\alias{do.pup.production}
\title{Pup production estimation...}
\usage{do.pup.production(fdir="")}

\description{
  Pup production estimation For Zc and Cu, creates data
  tables in the ACCESS database with the early pup
  mortality estimates during the season in each year.  Also
  constructs pdf plots of those values.
}

\details{
  Uses the \code{\link{production.stats}} to create
  ZcProduction and CuProduction tables for SMI in
  CIPinnipedCensusQuery.mdb. Any expected warnings or error
  messages are printed in PupProduction.log and unexpected
  errors/coding problems would be found in
  PupProduction.out.  If there are any errors then the
  production tables will not be created. In addition to
  creating the tables, it also produces a set of plots of
  the results in pdfs.
}
\value{None}
\author{Jeff Laake}
\arguments{
  \item{fdir}{directory for CIPinnipedCensusQuery.mdb}
}
