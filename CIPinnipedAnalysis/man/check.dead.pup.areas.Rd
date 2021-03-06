% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.dead.pup.areas.R
\name{check.dead.pup.areas}
\alias{check.dead.pup.areas}
\title{Checks for errors in the assigned dead pup sample areas and CU Survey dates}
\usage{
check.dead.pup.areas(fdir = NULL)
}
\arguments{
\item{fdir}{directory for CIPinnipedCensusMaster.mdb file}
}
\value{
None
}
\description{
Checks to make sure that each table in CIPinnipedCensusMaster.mdb has a
matching record in DeadPupSampleAreas table for SMI for species=Zc and Cu.
}
\details{
Also reports any CU survey dates in dead file that are more than 5 days from
assigned CU Survey Date for that survey number in that year.  Typically the error is

\preformatted{ If there are no missing records the output should look as
follows:
***Checking Zc CU dead pup census for Zc at SMI ***Checking Zc CU dead pup
census for Cu at SMI ***Checking Zc CU live pup census for Zc at SMI
***Checking Zc CU live pup census for Cu at SMI ***Checking Zc dead tag
initial for Zc at SMI ***Checking Zc dead tag resight for Zc at SMI } If
there are records that don't have a matching record in DeadPupSampleAreas
they are printed out. DeadPupSampleAreas must contain a record for each Area
Code for each year that the Area Code was used in either the tables for dead
pups or live pups because that tells the code how to sum the live and dead
counts for computation of the mortality rate.
}
\note{
Creates log file MissingDeadPupArea.txt with any errors
}
\author{
Jeff Laake
}
