

#' Pup production estimation
#' For Zc and Cu, creates data tables in the ACCESS database with the early pup
#' mortality estimates during the season in each year.  Also constructs pdf
#' plots of those values.
#' 
#' Uses the \code{\link{production.stats}} to create ZcProduction and
#' CuProduction tables for SMI in CIPinnipedCensusQuery.mdb. Any expected
#' warnings or error messages are printed in PupProduction.log and unexpected
#' errors/coding problems would be found in PupProduction.out.  If there are
#' any errors then the production tables will not be created. In addition to
#' creating the tables, it also produces a set of plots of the results in pdfs.
#' 
#' @param fdir directory for CIPinnipedCensusQuery.mdb
#' @return None
#' @export
#' @author Jeff Laake
do.pup.production=function(fdir="")
{
#
#  Constructs Pup Production tables for Cu and Zc at SMI in CIPinnipedCensusQuery
#  using live and dead pup counts from CIPinnipedCensusMaster via links and queries
#  in CIPinnipedCensusQuery.  It also saves plots of production in pdf files.
#
sink("PupProduction.log")
if(fdir=="")
{
	fdir1=system.file(package="CIPinnipedAnalysis")
	fdir2=fdir1
} else
{
	fdir1=fdir
    fdir2=file.path(fdir,"Master")
}
#  Make connection to CIPinnipedCensusQuery.mdb; assumed to be on J: (Calcur/Databases)
fdir1=file.path(fdir1,"CIPinnipedCensusQuery.mdb")
connection1=odbcConnectAccess2007(fdir1)
#  Make connection to CIPinnipedCensusMaster.mdb; assumed to be on J: (Calcur/Databases/Master)
fdir2=file.path(fdir2,"CIPinnipedCensusMaster.mdb")
connection2=odbcConnectAccess2007(fdir2)
# Delete current production tables	
xx=sqlDrop(connection1,"ZcProduction",errors=FALSE)
xx=sqlDrop(connection1,"CuProduction",errors=FALSE)
# Compute production stats for smi and castle rock for Zc
smidat=production.stats(island="SMI",mainland=TRUE,species="Zc",connection=connection2)
crdat=production.stats(island="SMI",mainland=FALSE,species="Zc",connection=connection2)
# Compute production stats for San Nicolas
livepups=sqlFetch(connection2,"Zc Cu live pup census")
snidat=production.stats(island="SNI",mainland=FALSE,species="Zc",connection=connection2,years=sort(unique(livepups$Year[livepups$Island=="SNI"])))
rm(livepups)
xx=sqlSave(connection1,rbind(smidat,crdat,snidat),tablename="ZcProduction",append=FALSE,rownames=FALSE)
# Compute production stats for CU on SMI and Castle Rock
smidat=production.stats(island="SMI",mainland=TRUE,species="Cu",connection=connection2)
crdat=production.stats(island="SMI",mainland=FALSE,species="Cu",connection=connection2)
xx=sqlSave(connection1,rbind(smidat,crdat),tablename="CuProduction",append=FALSE,rownames=FALSE)
sink()
# Get production tables and create pdfs
ZcProduction=sqlFetch(connection1,"ZcProduction")
CuProduction=sqlFetch(connection1,"CuProduction")
maxyr=max(ZcProduction$Year)

pdf("ZcPupProduction.pdf")
par(mfrow=c(2,1))
plot(ZcProduction$Year[ZcProduction$Area=="Mainland"&ZcProduction$Island=="SMI"],ZcProduction$PupProduction[ZcProduction$Area=="Mainland"&ZcProduction$Island=="SMI"],xlab="Year",ylab="Pup production",xaxt="n",main="SMI California sea lions Mainland")
lines(ZcProduction$Year[ZcProduction$Area=="Mainland"&ZcProduction$Island=="SMI"],ZcProduction$PupProduction[ZcProduction$Area=="Mainland"&ZcProduction$Island=="SMI"])
axis(1,at=seq(1970,maxyr,2),labels= seq(1970,maxyr,2)) 

zccr.years=ZcProduction$Year[ZcProduction$Area=="CastleRock"]
zccr.prod=ZcProduction$PupProduction[ZcProduction$Area=="CastleRock"]
zcyears=1970:max(zccr.years)
zcmissing.years=zcyears[!zcyears%in%zccr.years]
zccr.prod=c(zccr.prod,rep(NA,length(zcmissing.years)))
zccr.years=c(zccr.years,zcmissing.years)
zccr.prod=zccr.prod[order(zccr.years)]
zccr.years=sort(zccr.years)
plot(zccr.years,zccr.prod,main="California sea lions Castle Rock",xlab="Year",ylab="Pup production",xaxt="n",xlim=c(1970,maxyr))
lines(zccr.years,zccr.prod)
axis(1,at=seq(1970,maxyr,2),labels= seq(1970,maxyr,2))
dev.off()

pdf("ZcPupPreCensusMortality.pdf")
plot(ZcProduction$Year[ZcProduction$Area=="Mainland"&ZcProduction$Island=="SMI"],ZcProduction$MortalityRateAtLiveCount[ZcProduction$Area=="Mainland"&ZcProduction$Island=="SMI"],xlab="Year",ylab="Mortality",xaxt="n",main="Zc SMI Mainland Pre-census Mortality")
lines(ZcProduction$Year[ZcProduction$Area=="Mainland"&ZcProduction$Island=="SMI"],ZcProduction$MortalityRateAtLiveCount[ZcProduction$Area=="Mainland"&ZcProduction$Island=="SMI"])
axis(1,at=seq(1970,maxyr,2),labels= seq(1970,maxyr,2))
graphics.off()

pdf("ZcPupLiveCount.pdf")
par(mfrow=c(2,1))
plot(ZcProduction$Year[ZcProduction$Area=="Mainland"&ZcProduction$Island=="SMI"],ZcProduction$TotalLiveCount[ZcProduction$Area=="Mainland"&ZcProduction$Island=="SMI"],xlab="Year",ylab="Live pup count",xaxt="n",main="Mainland SMI California sea lion Live Pup Count")
lines(ZcProduction$Year[ZcProduction$Area=="Mainland"&ZcProduction$Island=="SMI"],ZcProduction$TotalLiveCount[ZcProduction$Area=="Mainland"&ZcProduction$Island=="SMI"])
axis(1,at=seq(1970,maxyr,2),labels= seq(1970,maxyr,2))

zccr.years=ZcProduction$Year[ZcProduction$Area=="CastleRock"]
zccr.prod=ZcProduction$PupProduction[ZcProduction$Area=="CastleRock"]
zcyears=1970:max(zccr.years)
zcmissing.years=zcyears[!zcyears%in%zccr.years]
zccr.prod=c(zccr.prod,rep(NA,length(zcmissing.years)))
zccr.years=c(zccr.years,zcmissing.years)
zccr.prod=zccr.prod[order(zccr.years)]
zccr.years=sort(zccr.years)
plot(zccr.years,zccr.prod,main="Castle Rock California sea lion Live Pup Count",xlab="Year",ylab="Live pup count",xaxt="n",xlim=c(1970,maxyr))
lines(zccr.years,zccr.prod)
axis(1,at=seq(1970,maxyr,2),labels= seq(1970,maxyr,2))
dev.off()

pdf("CuPupProduction.pdf")
par(mfrow=c(2,1))
plot(CuProduction$Year[CuProduction$Area=="Mainland"],CuProduction$PupProduction[CuProduction$Area=="Mainland"],xlab="Year",ylab="Pup production",xaxt="n",main="Northern fur seals Mainland")
lines(CuProduction$Year[CuProduction$Area=="Mainland"],CuProduction$PupProduction[CuProduction$Area=="Mainland"])
axis(1,at=seq(1970,maxyr,2),labels= seq(1970,maxyr,2))


cr.years=CuProduction$Year[CuProduction$Area=="CastleRock"]
cr.prod=CuProduction$PupProduction[CuProduction$Area=="CastleRock"]
years=1970:max(cr.years)
missing.years=years[!years%in%cr.years]
cr.prod=c(cr.prod,rep(NA,length(missing.years)))
cr.years=c(cr.years,missing.years)
cr.prod=cr.prod[order(cr.years)]
cr.years=sort(cr.years)
plot(cr.years,cr.prod,main="Northern fur seals Castle Rock",xlab="Year",ylab="Pup production",xaxt="n",xlim=c(1970,maxyr))
lines(cr.years,cr.prod)
axis(1,at=seq(1970,maxyr,2),labels= seq(1970,maxyr,2))
dev.off()

pdf("CuPupPreCensusMortality.pdf")
plot(CuProduction$Year[CuProduction$Area=="Mainland"],CuProduction$MortalityRateAtLiveCount[CuProduction$Area=="Mainland"],xlab="Year",ylab="Mortality",xaxt="n",main="Cu Mainland Pre-census Mortality",xlim=c(1970,maxyr))
lines(CuProduction$Year[CuProduction$Area=="Mainland"],CuProduction$MortalityRateAtLiveCount[CuProduction$Area=="Mainland"])
axis(1,at=seq(1970,maxyr,2),labels= seq(1970,maxyr,2))
dev.off()

pdf("CuPupLiveCount.pdf")
par(mfrow=c(2,1))
plot(CuProduction$Year[CuProduction$Area=="Mainland"],CuProduction$TotalLiveCount[CuProduction$Area=="Mainland"],xlab="Year",ylab="Live pup count",xaxt="n",main="Mainland Northern fur seal Live Pup Count",xlim=c(1970,maxyr),ylim=c(0,2500))
lines(CuProduction$Year[CuProduction$Area=="Mainland"],CuProduction$TotalLiveCount[CuProduction$Area=="Mainland"])
axis(1,at=seq(1970,maxyr,2),labels= seq(1970,maxyr,2))

cr.years=CuProduction$Year[CuProduction$Area=="CastleRock"]
cr.prod=CuProduction$PupProduction[CuProduction$Area=="CastleRock"]
years=1970:max(cr.years)
missing.years=years[!years%in%cr.years]
cr.prod=c(cr.prod,rep(NA,length(missing.years)))
cr.years=c(cr.years,missing.years)
cr.prod=cr.prod[order(cr.years)]
cr.years=sort(cr.years)
plot(cr.years,cr.prod,main="Castle Rock Northern fur seal Live Pup Count",xlab="Year",ylab="Pup production",xaxt="n",xlim=c(1970,maxyr), ylim=c(0,1500))
lines(cr.years,cr.prod)
axis(1,at=seq(1970,maxyr,2),labels= seq(1970,maxyr,2))
dev.off()

odbcCloseAll()
invisible()
}

