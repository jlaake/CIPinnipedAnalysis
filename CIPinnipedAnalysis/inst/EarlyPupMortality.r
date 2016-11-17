#' Early pup mortality estimation
#' 
#' For Zc and Cu, creates data tables in the ACCESS database with the early pup
#' mortality estimates during the season in each year.  Also constructs pdf
#' plots of those values.
#' 
#' Uses the \code{\link{mortality.stats}}function to create two data tables of
#' early pup mortality for Cu and Zc in CIPinnipedCensusQuery.mdb and it
#' creates plots in pdf files. If there are any errors then the mortality
#' tables will not be created. If successful the pdf files are created in this
#' directory which can then be copied to Presentations/PinnipedPlots.
#' 
#'
#'  Produces early pup mortality estimates for Cu and Zc at SMI using
#'  live and dead counts in CIPinnipedCensus database.
#'  A table is constructed for each species in CIPinnipedCensusQuery.mdb
#'  and plots of the cummulative survival by year are created and save as
#'  pdf files
#
if(!exists("fdir"))fdir=NULL
#
#  Setup directories for locations of data
#
if(!is.null(fdir) && fdir=="")
{
	fdir1=system.file(package="CalcurData")
	fdir2=fdir1
} else
{
	fdir1=fdir
	fdir2=fdir
}
sdir=system.file(package="CIPinnipedAnalysis")
source(file.path(sdir,"ZcEarlyPupMortality.r"))
source(file.path(sdir,"CuEarlyPupMortality.r"))
