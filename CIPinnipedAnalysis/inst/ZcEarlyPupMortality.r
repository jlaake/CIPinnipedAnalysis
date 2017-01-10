' Early pup mortality estimation
#' 
#' For Zc creates data tables in the ACCESS database with the early pup
#' mortality estimates during the season in each year.  Also constructs pdf
#' plots of those values.
#' 
#' Uses the \code{\link{mortality.stats}}function to create data table of
#' early pup mortality for  Zc in CIPinnipedCensusQuery.mdb and it
#' creates plots in pdf files. If there are any errors then the mortality
#' tables will not be created. If successful the pdf files are created in this
#' directory which can then be copied to Presentations/PinnipedPlots.
#' 
#'
#'  Produces early pup mortality estimates for  Zc at SMI using
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
#Set maxyear for SMI
deads=getCalcurData("CIPCensus","Zc Cu dead pup census")
maxyear=max(deads$Year[deads$Island=="SMI"])

# Construct mortality tables for Zc on San Miguel
zcsmi.mort=NULL
Production=getCalcurData("CIPquery","ZcProduction",dir=fdir1)
Production=Production[Production$Area=="Mainland"&Production$Island=="SMI",]
zcsmi.mort=NULL
for(y in 1991:maxyear)
{
	dead=correct_dead("SMI",y)$bystrata
	pups=Production$AdjustedDeadInDeadSampleArea[Production$Year==y]+Production$LiveInDeadSampleArea[Production$Year==y]
	zcsmi.mort=rbind(zcsmi.mort,MortalityStats(year=y,dead,pups))
}	
zcsmi.mort$Island="SMI"


#Set maxyear for SNI
maxyear=max(deads$Year[deads$Island=="SNI"])
Production=getCalcurData("CIPquery","ZcProduction",dir=fdir1)
Production=Production[Production$Island=="SNI",]
zcsni.mort=NULL
for(y in 2005:maxyear)
{
	dead=correct_dead("SNI",y)$bystrata
	pups=Production$AdjustedDeadInDeadSampleArea[Production$Year==y]+Production$LiveInDeadSampleArea[Production$Year==y]
	if(y<2008)pups=NA
	zcsni.mort=rbind(zcsni.mort,MortalityStats(year=y,dead,pups))
}	
zcsni.mort$Island="SNI"
xx=saveCalcurData(rbind(zcsmi.mort,zcsni.mort),db="CIPquery",tbl="ZcEarlyPupMortality",dir=fdir1)

