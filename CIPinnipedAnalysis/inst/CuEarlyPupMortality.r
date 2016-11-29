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


# Construct mortality tables for Cu on San Miguel
cu.mort=CuMortalityStats(fdir1=fdir1,fdir2=fdir2)
cu.mort.table=cu.mort
cu.mort.table$SurveyDate=substr(as.character(cu.mort$SurveyDate),1,10)
xx=saveCalcurData(cu.mort.table,db="CIPquery",tbl="CuEarlyPupMortality",dir=fdir1)
pdf("CuEarlyPupMortality.pdf",width=9)
minyear=min(cu.mort$Year)
maxyear=max(cu.mort$Year)
numyears=maxyear-minyear+1
ltype=rep(1:6,ceiling(numyears/6))
col=rep(c("black","red","blue","orange","green"),each=ceiling(numyears/5))
cu.mort$DaysFrom15June=as.Date(cu.mort$SurveyDate)-as.Date(paste(cu.mort$Year,"-06-15",sep=""))
i=1
with(cu.mort[cu.mort$Year==minyear,],plot(c(0,DaysFrom15June),c(1,CumS),xlab="Days from 15 June",ylab="Cumulative Survival",ylim=c(min(cu.mort$CumS,na.rm=TRUE),1),lty=ltype[1],col=col[1],type="l",xlim=c(0,max(cu.mort$DaysFrom15June))))
for (y in (minyear+1):maxyear)
{
	i=i+1
	with(cu.mort[cu.mort$Year==y,],lines(c(0,DaysFrom15June),c(1,CumS),lty=ltype[i],col=col[i]))
}
legend(140,1,legend=minyear:maxyear,lty=ltype[1:numyears],col=col[1:numyears],cex=0.8)

CumSto1Aug=rep(NA,numyears)
i=1
for (y in minyear:maxyear)
{
	cat("\n",y)
	CumSto1Aug[i]=with(cu.mort[cu.mort$Year==y,],approx(DaysFrom15June,CumS,xout=42)$y)
	i=i+1
}
plot(minyear:maxyear,CumSto1Aug,xlab="Year",ylab="Cumulative survival to 27 July",type="b")

dev.off()



