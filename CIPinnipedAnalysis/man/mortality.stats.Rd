\name{mortality.stats}
\alias{mortality.stats}
\title{Performs numerical computations for early pup mortality estimation}
\usage{
  mortality.stats(PreLiveCountCf = 1.33,
    PostLiveCountCf = 1.25, species = "Zc", island = "SMI",
    fdir1, fdir2, years = NULL)
}
\arguments{
  \item{PreLiveCountCf}{Multiplicative correction factor
  for observed mortality prior to live count to account for
  dead pups that were missed and decomposed or were buried}

  \item{PostLiveCountCf}{Multiplicative correction factor
  for observed mortality after live count}

  \item{species}{either "Zc" for Zalophus or "Cu" for
  Callorhinus}

  \item{island}{the character name for the island "SMI" or
  "SNI"}

  \item{fdir1}{database directory for
  CIPinnipedCensusQuery}

  \item{fdir2}{database directory for
  CIPinnipedCensusMaster}

  \item{years}{vector of years that should be selected if
  not all years}
}
\value{
  dataframe containing year,SurveyDate,Pups(number of pups
  alive at time),ObservedDead (number observed to die in
  that interval) AdjustedDead(corrected count that
  died),MortalityRate(mortality rate during the interval
  between surveys),DailyMortalityRate (mortality rate on a
  constant daily rate), CumS (cummulative survival to that
  time).
}
\description{
  Creates a survivorship curve from dead pup count surveys.
  The curve extends from birth to either Aug or Sept-Oct
  depending on the extent of the mortality surveys which
  has varied across the years.
}
\details{
  \code{\link{do.early.pup.mortality}} is the driver
  function to construct the ZcEarlPupMortality and
  CuEarlyPupMortality tables in CIPinnipedCensusQuery.mdb.
  This function uses the ZcProduction and CuProduction
  tables created by \code{\link{do.pup.production}} and
  should be run after it is updated with any new data.
}
\author{
  Jeff Laake
}

