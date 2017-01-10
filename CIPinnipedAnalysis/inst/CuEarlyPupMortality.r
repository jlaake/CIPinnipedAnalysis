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

#Set maxyear for SMI
deads=getCalcurData("CIPCensus","Zc Cu dead pup census")
maxyear=max(deads$Year[deads$Island=="SMI"])

# Construct mortality tables for Cu on San Miguel
Production=getCalcurData("CIPquery","CuProduction",dir=fdir1)
Production=Production[Production$Area=="Mainland"&Production$Island=="SMI",]
cusmi.mort=NULL
for(y in 1997:maxyear)
{
	dead=cu_correct_dead(y)$bystrata
	pups=Production$AdjustedDeadInDeadSampleArea[Production$Year==y]+Production$LiveInDeadSampleArea[Production$Year==y]
	cusmi.mort=rbind(cusmi.mort,MortalityStats(year=y,dead,pups))
}	
cusmi.mort$Island="SMI"
cu.mort.table=cusmi.mort
xx=saveCalcurData(cu.mort.table,db="CIPquery",tbl="CuEarlyPupMortality",dir=fdir1)





