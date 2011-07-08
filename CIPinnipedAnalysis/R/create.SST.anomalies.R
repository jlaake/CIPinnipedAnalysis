

#' Create SST monthly anomaly matrices for each of 8 locations
#' Computes SST anomaly matrix for month and year for each of 8 locations: East
#' Santa Barbara (ESB), West Santa Barbara (WSB), Point Arguello (PtArg), Point
#' Santa Maria (PtSM), Port San Luis (PtSL), Cape San Martin (CSM), Monterey
#' Bay (MB) and Point Reyes (PtReyes).
#' 
#' Creates daily means from hourly observations.  Averages daily means within
#' month.  Computes average monthly value from a defined set of years
#' (\code{average.years}) and subtracts the monthly means from the observed
#' year-month matrix to produce anomaly matrix.
#' 
#' @param average.years vector of years to be included in the monthly mean
#'   calculations the monthly means are subtracted from the specific year-month
#'   values to create the anomaly
#' @param fdir directory for location of environmental.mdb file. default
#'   assumes it is in directory containing the workspace
#' @param store if TRUE, deletes and replaces anaomaly tables in the ACCESS
#'   database
#' @return A three dimensional array (year by month by location) with the 8
#'   locations ordered as ESB, WSB, PtArg, PtSM, PtSL, CSM, MB and PtReyes.
#'   Side effect of function is to delete existing anomaly tables in ACCESS and
#'   re-write new ones.
#' @export
#' @author Jeff Laake
#' @seealso \code{\link{create.anomalies}}
#' @examples
#' 
#' # fdir should be set to location of environmental.data.mdb
#' # each defaults below to current working directory
#' # Use "J:/" if it should be Calcurr/Databases
#' # fdir="J:/"
#' fdir="" 
#' edir=""
#' anomalies=create.SST.anomalies(c(1994:1996,1998:2008),fdir=fdir)
#' CentralSSTAnomalies=t(apply(anomalies[,,3:7],c(2,1),mean,na.rm=TRUE))
#' SouthSSTAnomalies=t(apply(anomalies[,,1:2],c(2,1),mean,na.rm=TRUE))
#' NorthSSTAnomalies=anomalies[,,8]
#' SSTAnomalies=t(apply(anomalies[,,1:5],c(2,1),mean,na.rm=TRUE))
#' maxyear= max(as.numeric(row.names(SSTAnomalies)))
#' minyear= min(as.numeric(row.names(SSTAnomalies)))
#' numyears=maxyear-minyear+1
#' JantoMayAnomalies=rowMeans(SSTAnomalies[,c("Jan","Feb","Mar","Apr","May")])
#' JunetoSeptAnomalies=rowMeans(SSTAnomalies[,c("June","July","Aug","Sept")])
#' OcttoDecAnomalies=rowMeans(SSTAnomalies[,c("Oct","Nov","Dec")])
#' x=rbind(data.frame(Year=minyear:maxyear,Season=rep("Spring",numyears),SSTAnomaly=JantoMayAnomalies),
#'   data.frame(Year=minyear:maxyear,Season=rep("Summer",numyears),SSTAnomaly=JunetoSeptAnomalies),
#'   data.frame(Year=minyear:maxyear,Season=rep("Fall",numyears),SSTAnomaly=OcttoDecAnomalies) )
#' x$Season=factor(x$Season,levels=c("Spring","Summer","Fall"))
#' x=x[order(x$Year,x$Season),]
#' pdf("SSTAnomaly.pdf",pointsize=10)
#' par(mfrow=c(3,1))
#' plot(minyear:maxyear,x$SSTAnomaly[x$Season=="Spring"],ylab="SST Anomaly",main="Winter/Spring (Jan-May)",xlab="")
#' lines(minyear:maxyear,x$SSTAnomaly[x$Season=="Spring"])
#' abline(h=0)
#' plot(minyear:maxyear,x$SSTAnomaly[x$Season=="Summer"],ylab="SST Anomaly",main="Summer (June-Sept)",xlab="")
#' lines(minyear:maxyear,x$SSTAnomaly[x$Season=="Summer"])
#' abline(h=0)
#' plot(minyear:maxyear,x$SSTAnomaly[x$Season=="Fall"],ylab="SST Anomaly",main="Fall (Oct-Dec)",xlab="")
#' lines(minyear:maxyear,x$SSTAnomaly[x$Season=="Fall"])
#' abline(h=0)
#' dev.off()
#' 
#' fpath=file.path(system.file(package="CIPinnipedAnalysis"),"environmental.data.mdb")
#' require(RODBC)
#' connection=odbcConnectAccess2007(fpath)
#' pdf("MultivariateENSOIndex.pdf",pointsize=10)
#' MEI=sqlFetch(connection,"MEI")
#' minyear=min(MEI$Year)
#' maxyear=max(MEI$Year)
#' numyears=maxyear-minyear+1
#' plot(MEI$MEI,type="l",lwd=2,xaxt="n",ylab="MEI",xlab="Year")
#' axis(1,at=12*(0:(numyears-1))+1,labels=as.character(minyear:maxyear))
#' abline(h=0)
#' abline(h=1)
#' abline(h=-1)
#' dev.off()
#' 
#' pdf("UWIAnomaly.pdf",pointsize=10)
#' par(mfrow=c(2,1))
#' UWI=sqlFetch(connection,"UWIAnomaly")
#' UWI=UWI[order(UWI$Year,UWI$Month),]
#' minyear=min(UWI$Year)
#' maxyear=max(UWI$Year)
#' numyears=maxyear-minyear+1
#' plot(UWI$UWI[UWI$Location=="36N122W"],type="l",lwd=2,xaxt="n",ylab="UWI",xlab="Year",main="36N122W")
#' axis(1,at=12*(0:(numyears-1))+1,labels=as.character(minyear:maxyear))
#' abline(h=0)
#' plot(UWI$UWI[UWI$Location=="33N119W"],type="l",lwd=2,xaxt="n",ylab="UWI",xlab="Year",main="33N119W")
#' axis(1,at=12*(0:(numyears-1))+1,labels=as.character(minyear:maxyear))
#' abline(h=0)
#' dev.off()
#' 
#' odbcClose(connection)
#' 
create.SST.anomalies=function(average.years,fdir="",store=FALSE)
{
# Create SST monthly anomaly matrices by location
#
# Arguments:
#
#  average.years - vector of years to be included in the monthly mean calculations
#                   the monthly means are subtracted from the specific year-month values to create the anomaly
#  dir           - directory for location of environmental.mdb file. default assumes it is in directory containing the
#                   workspace
#  store         - if TRUE, stores new anomaly datafiles in ACCESS
#
# Value:
#   returns a three dimensional array (year by month by location) with the 8 locations
#    ordered as ESB, WSB, PtArg, PtSM, PtSL, CSM, MB and PtReyes.
#
# Side effect of function is to delete existing anomaly tables in ACCESS and re-write new ones.
#
# Open a connection to the ACCESS database, attach the data tables and
# modify some fieldnames and values to make uniform
  if(fdir=="")fdir=system.file(package="CIPinnipedAnalysis")
  fdir=file.path(fdir,"environmental.data.mdb")
  connection=odbcConnectAccess2007(fdir)
  PtArg=sqlFetch(connection,"PtArguelloBuoyData")
  PtSm=sqlFetch(connection,"PtSantaMariaBuoyData")
  ESB=sqlFetch(connection,"EastSantaBarbaraChannelBuoyData")
  WSB=sqlFetch(connection,"WestSantaBarbaraChannelBuoyData")
  PtSanLuis=sqlFetch(connection,"PtSanLuisDailySST")
  Nsst=sqlFetch(connection,"NorthernCaliforniaSSTData")
  Nsst$YYYY=as.POSIXlt(Nsst$Date)$year+1900
  Nsst$date=Nsst$Date
  Nsst$WTMP=Nsst$SST
  names(PtSanLuis)=c("date","WTMP")
# Compute full range of years for computation of anomaly matrices
  years=c(PtArg$YYYY,PtSm$YYYY,ESB$YYYY,WSB$YYYY,as.POSIXlt(PtSanLuis$date)$year+1900,Nsst$YYYY)
  years=years[years>1900]
  yrange=min(years,na.rm=TRUE):max(years,na.rm=TRUE)
# Compute anomaly matrix for each location
  PtArgAnomalies=create.anomalies(PtArg,yrange,average.years)
  anomalies.bylocation=array(NA,dim=c(dim(PtArgAnomalies),8))
  anomalies.bylocation[,,1]=create.anomalies(ESB,yrange,average.years)
  anomalies.bylocation[,,2]=create.anomalies(WSB,yrange,average.years)
  anomalies.bylocation[,,3]=PtArgAnomalies
  anomalies.bylocation[,,4]=create.anomalies(PtSm,yrange,average.years)
  anomalies.bylocation[,,5]=create.anomalies(PtSanLuis,yrange,average.years)
  anomalies.bylocation[,,6]=create.anomalies(Nsst[Nsst$Location=="Cape San Martin",],yrange,average.years)
  anomalies.bylocation[,,7]=create.anomalies(Nsst[Nsst$Location=="Monterey Bay",],yrange,average.years)
  anomalies.bylocation[,,8]=create.anomalies(Nsst[Nsst$Location=="Pt Reyes",],yrange,average.years)
  dimnames(anomalies.bylocation)=list(yrange,c("Jan","Feb","Mar","Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"),
                                        c("ESB","WSB","PtArg","PtSM","PtSL","CSM","MB","PtReyes"))
  if(store)
  {
# Compute anomaly dataframes and save in ACCESS
     PtArguelloSSTAnomalies=create.anomaly.table(PtArgAnomalies)
	 PtSanLuisSSTAnomalies=create.anomaly.table( anomalies.bylocation[,,5])
	 PtSantaMariaSSTAnomalies=create.anomaly.table( anomalies.bylocation[,,4])
	 WestSantaBarbaraSSTAnomalies=create.anomaly.table( anomalies.bylocation[,,2])
	 EastSantaBarbaraSSTAnomalies=create.anomaly.table( anomalies.bylocation[,,1])
	 MontereyBaySSTAnomalies=create.anomaly.table( anomalies.bylocation[,,7])
	 CSMSSTAnomalies=create.anomaly.table( anomalies.bylocation[,,6])
	 PtReyesSSTAnomalies=create.anomaly.table( anomalies.bylocation[,,8])
	 xx=sqlDrop(connection,"PtSanLuisSSTAnomalies",errors=FALSE)
	 xx=sqlDrop(connection,"PtArguelloSSTAnomalies",errors=FALSE)
	 xx=sqlDrop(connection,"PtSantaMariaSSTAnomalies",errors=FALSE)
	 xx=sqlDrop(connection,"WestSantaBarbaraSSTAnomalies",errors=FALSE)
	 xx=sqlDrop(connection,"EastSantaBarbaraSSTAnomalies",errors=FALSE)
	 xx=sqlDrop(connection,"MontereyBaySSTAnomalies",errors=FALSE)
	 xx=sqlDrop(connection,"CapeSanMartinSSTAnomalies",errors=FALSE)
	 xx=sqlDrop(connection,"PtReyesSSTAnomalies",errors=FALSE)
	 xx=sqlSave(connection,PtArguelloSSTAnomalies,tablename="PtArguelloSSTAnomalies",append=FALSE,rownames=FALSE)
	 xx=sqlSave(connection,PtSanLuisSSTAnomalies,tablename="PtSanLuisSSTAnomalies",append=FALSE,rownames=FALSE)
	 xx=sqlSave(connection,PtSantaMariaSSTAnomalies,tablename="PtSantaMariaSSTAnomalies",append=FALSE,rownames=FALSE)
	 xx=sqlSave(connection,WestSantaBarbaraSSTAnomalies,tablename="WestSantaBarbaraSSTAnomalies",append=FALSE,rownames=FALSE)
	 xx=sqlSave(connection,EastSantaBarbaraSSTAnomalies,tablename="EastSantaBarbaraSSTAnomalies",append=FALSE,rownames=FALSE)
	 xx=sqlSave(connection,MontereyBaySSTAnomalies,tablename="MontereyBaySSTAnomalies",append=FALSE,rownames=FALSE)
	 xx=sqlSave(connection,CSMSSTAnomalies,tablename="CapeSanMartinSSTAnomalies",append=FALSE,rownames=FALSE)
	 xx=sqlSave(connection,PtReyesSSTAnomalies,tablename="PtReyesSSTAnomalies",append=FALSE,rownames=FALSE)
  }
  odbcClose(connection)
  return(anomalies.bylocation)
}
