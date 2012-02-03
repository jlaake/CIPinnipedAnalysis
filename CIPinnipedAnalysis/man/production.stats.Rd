\name{production.stats}
\alias{production.stats}
\title{Pup production stats
: computes pup production stats for Zalophus or Callorhinus at San Miguel or
Castle Rock and creates production tables in the ACCESS database.}
\usage{
  production.stats(island = "SMI", mainland = TRUE,
    years = NULL, PreLiveCountCf = 1.33, species = "Zc",
    connection)
}
\arguments{
  \item{island}{the character name for the island "SMI" or
  "SNI"}

  \item{mainland}{If TRUE, computation is for San Miguel
  mainland; otherwise, Castle Rock}

  \item{years}{vector of years to select or NULL if all to
  be used}

  \item{PreLiveCountCf}{Multiplicative correction factor
  for observed mortality prior to live count to account for
  dead pups that were missed and decomposed or were buried}

  \item{species}{either "Zc" for Zalophus or "Cu" for
  Callorhinus}

  \item{connection}{database connection in CensusMaster}
}
\value{
  dataframe containing results
  LiveCountDate,Area,Year,LiveInDeadSampleArea,DeadInDeadSampleArea,AdjustedDeadInDeadSampleArea,
  MortalityRateAtLiveCount,TotalLiveCountByYear,PupProduction
}
\description{
  Dead pup surveys are conducted typically on 3 occasions
  during July.  Dead pups are either stacked or tagged so
  they are not recounted.  A live count is also conducted
  in late July and it often occurs between two of the dead
  pup surveys.  This function totals up the dead pups,
  applies an average multiplicative correction factor to
  account for pups that died and were never counted, and
  then interpolates the number that died up to the time of
  the live count.  The dead pup surveys are only conducted
  in a portion of the rookery at SMI so the estimated
  number of dead pups in the dead pup sample area and the
  live count in that same area are used to create a
  pre-census mortality rate.  That rate is then used to
  scale the total live count to compute the number of pups
  produced (born). The computed values are output in the
  production table for each species.
}
\details{
  \code{\link{do.pup.production}} is the driver function
  which calls this function for Zc for SMI and Castle Rock
  and then combines the result in the table ZcProduction
  which is stored in CIPinnipedCensusQuery.mdb.  Then it
  does the same thing for Cu and creates CuProduction.
}
\author{
  Jeff Laake
}

