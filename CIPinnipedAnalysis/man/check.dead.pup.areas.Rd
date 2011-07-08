\name{check.dead.pup.areas}
\alias{check.dead.pup.areas}
\title{Checks for errors in the assigned dead pup sample areas...}
\usage{check.dead.pup.areas(fdir="")}

\description{
  Checks for errors in the assigned dead pup sample areas
  Checks to make sure that each table in
  CIPinnipedCensusMaster.mdb has a matching record in
  DeadPupSampleAreas table for SMI for species=Zc and Cu.
}

\details{
  \preformatted{ If there are no missing records the output
  should look as follows:

  ***Checking Zc CU dead pup census for Zc at SMI
  ***Checking Zc CU dead pup census for Cu at SMI
  ***Checking Zc CU live pup census for Zc at SMI
  ***Checking Zc CU live pup census for Cu at SMI
  ***Checking Zc dead tag initial for Zc at SMI ***Checking
  Zc dead tag resight for Zc at SMI } If there are records
  that don't have a matching record in DeadPupSampleAreas
  they are printed out. DeadPupSampleAreas must contain a
  record for each Area Code for each year that the Area
  Code was used in either the tables for dead pups or live
  pups because that tells the code how to sum the live and
  dead counts for computation of the mortality rate.
}
\value{None}
\note{Creates log file MissingDeadPupArea.txt with any errors}
\author{Jeff Laake}
\arguments{
  \item{fdir}{directory for CIPinnipedCensusMaster.mdb file}
}
