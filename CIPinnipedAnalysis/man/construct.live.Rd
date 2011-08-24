\name{construct.live}
\alias{construct.dead}
\alias{construct.live}
\title{Constructs dataframes of live or dead pup counts
:extracts data from CIPinnipedMaster database and constructs dataframe of
dead pup counts or live pup counts for specified island and species.}
\usage{
  construct.live(island, species, connection)
}
\arguments{
  \item{island}{Either "SMI" or "SNI"}

  \item{species}{Either "Zc" or "Cu"}

  \item{connection}{RODBC connection to the database}
}
\value{
  for \code{construct.live} a dataframe containing
  YearArea, Survey date, Area, AvgCount, Year, Month, Day,
  DeadPupArea (Y/N if included in area sampled for dead
  pups).

  for \code{construct.dead} a dataframe containing
  YearArea, SurveyNumber, SurveyDate, Area, count, Year,
  Month, Day, DeadPupArea (Y/N if included in area sampled
  for dead pups), and MatchingLiveArea.
}
\description{
  \code{construct.live} extracts data from table \code{Zc
  Cu live pup census} and sums counts for an observer
  within an area for a year and then averages those counts
  across observers that counted the same area that year.
}
\details{
  \code{construct.dead} extracts data from table \code{Zc
  Cu dead pup census} and table \code{Zc dead tag initial}
  and totals dead pups by survey number/date.  It will
  produce errors if there are more than one survey date for
  each survey number within an area in a year.
}
\author{
  Jeff Laake
}

