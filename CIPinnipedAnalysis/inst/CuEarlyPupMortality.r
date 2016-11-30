' Early pup mortality estimation
#' 
#' For  Cu, creates data tables in the ACCESS database with the early pup
#' mortality estimates during the season in each year.  Also constructs pdf
#' plots of those values.
#' 
#' Uses the \code{\link{mortality.stats}}function to creates data table of
#' early pup mortality for Cu in CIPinnipedCensusQuery.mdb and it
#' creates plots in pdf files. If there are any errors then the mortality
#' tables will not be created. If successful the pdf files are created in this
#' directory which can then be copied to Presentations/PinnipedPlots.
#' 
#'
#'  Produces early pup mortality estimates for Cu at SMI using
#'  live and dead counts in CIPinnipedCensus database.
#'  A table is constructed for each species in CIPinnipedCensusQuery.mdb
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


# Construct mortality tables for Cu on San Miguel
cu.mort=CuMortalityStats(fdir1=fdir1,fdir2=fdir2)
cu.mort.table=cu.mort
cu.mort.table$SurveyDate=substr(as.character(cu.mort$SurveyDate),1,10)
xx=saveCalcurData(cu.mort.table,db="CIPquery",tbl="CuEarlyPupMortality",dir=fdir1)


