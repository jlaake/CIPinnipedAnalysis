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
#' @param fdir directory for data files; if NULL uses location specified in databases.txt of CalcurData package; if "" uses databases in CalcurData pacakge; otherwise uses specified directory location
#' @return None
#' @export
#' @author Jeff Laake
do.early.pup.mortality=function(fdir=NULL)
{
#
#  Produces early pup mortality estimates for Cu and Zc at SMI using
#  live and dead counts in CIPinnipedCensus database.
#  A table is constructed for each species in CIPinnipedCensusQuery.mdb
#  and plots of the cummulative survival by year are created and save as
#  pdf files
#
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
# Construct mortality tables for Zc on San Miguel
zcsmi.mort=mortality.stats(island="SMI",species="Zc",fdir1=fdir1,fdir2=fdir2)
zcsmi.mort$Island="SMI"
#zcsni.mort=mortality.stats(island="SNI",species="Zc",fdir1=fdir1,fdir2=fdir2,years=2005:2006)
#zcsni.mort$Island="SNI"
#zc.mort=rbind(zcsmi.mort,zcsni.mort)
zc.mort=zcsmi.mort
# Construct mortality tables for Cu on San Miguel
cu.mort=mortality.stats(species="Cu",fdir1=fdir1,fdir2=fdir2)
cu.mort.table=cu.mort
zc.mort.table=zc.mort
cu.mort.table$SurveyDate=substr(as.character(cu.mort$SurveyDate),1,10)
zc.mort.table$SurveyDate=substr(as.character(zc.mort$SurveyDate),1,10)
xx=saveCalcurData(cu.mort.table,db="CIPquery",tbl="CuEarlyPupMortality",dir=fdir1)
xx=saveCalcurData(zc.mort.table,db="CIPquery",tbl="ZcEarlyPupMortality",dir=fdir1)
pdf("CuEarlyPupMortality.pdf",width=9)
minyear=min(cu.mort$Year)-1
maxyear=max(cu.mort$Year)
numyears=maxyear-minyear
for (i in 1:ceiling(numyears/6))
{
  par(mfrow=c(2,3))
  for( year in (minyear+(i-1)*6+1):min(maxyear,minyear+i*6))
  {
     x=c(as.POSIXct(paste(year,"-06-15",sep="")),cu.mort$SurveyDate[cu.mort$Year==year])
     xlim=c(as.POSIXct(paste(year,"-06-15",sep="")),as.POSIXct(paste(year,"-10-15",sep="")))
     y=c(1,cu.mort$CumS[cu.mort$Year==year])
     plot(x,y,main=year,xlab="Date",ylab="Cumulative Survival",ylim=c(min(cu.mort$CumS),1),type="b",xlim=xlim)
  }
}
dev.off()

pdf("ZcEarlyPupMortality.pdf",width=9)
zc.mort=zc.mort[zc.mort$Island=="SMI",]
minyear=min(zc.mort$Year)-1
maxyear=max(zc.mort$Year)
numyears=maxyear-minyear
for (i in 1:ceiling(numyears/6))
{
  par(mfrow=c(2,3))
  for( year in (minyear+(i-1)*6+1):min(maxyear,minyear+i*6))
  {
     x=c(as.POSIXct(paste(year,"-06-15",sep="")),zc.mort$SurveyDate[zc.mort$Year==year])
     y=c(1,zc.mort$CumS[zc.mort$Year==year])
     xlim=c(as.POSIXct(paste(year,"-06-15",sep="")),as.POSIXct(paste(year,"-10-15",sep="")))
     plot(x,y,main=year,xlab="Date",ylab="Cumulative Survival",ylim=c(min(zc.mort$CumS),1),type="b",xlim=xlim)
     }
}
dev.off()
invisible()
}

