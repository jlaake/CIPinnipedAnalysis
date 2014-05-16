\name{plot.weight.series}
\alias{plot.weight.series}
\title{Plot weight series with confidence interval}
\usage{
\method{plot}{weight.series}(time, predictions, date, ...)
}
\arguments{
  \item{time}{vector of years}

  \item{predictions}{list with elements fit (estimate) and
  se (std error)}

  \item{date}{string for title for date weights are
  predicted}

  \item{...}{other arguments for plot}
}
\value{
None
}
\description{
Plots time series of predicted weights and conf intervals
}
\author{
Jeff Laake
}

