#' Create SST monthly anomaly matrices for each of 8 locations
#' 
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
#' @param fdir directory for data files; if NULL uses location specified in databases.txt of CalcurData package; if "" uses databases in CalcurData pacakge; otherwise uses specified directory location 
#' @param store if TRUE, deletes and replaces anaomaly tables in the ACCESS
#'   database
#' @return A three dimensional array (year by month by location) with the 8
#'   locations ordered as ESB, WSB, PtArg, PtSM, PtSL, CSM, MB and PtReyes.
#'   Side effect of function is to delete existing anomaly tables in ACCESS and
#'   re-write new ones.
#' @export
#' @author Jeff Laake
#' @seealso \code{\link{create.anomalies}}
#' 
create.SST.anomalies=function(average.years,fdir=NULL,store=FALSE)
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
    if(fdir=="")fdir=system.file(package="CalcurData")
    PtArg=getCalcurData("Environ","PtArguelloBuoyData",dir=fdir)
	PtSm=getCalcurData("Environ","PtSantaMariaBuoyData",dir=fdir)
	ESB=getCalcurData("Environ","EastSantaBarbaraChannelBuoyData",dir=fdir)
	WSB=getCalcurData("Environ","WestSantaBarbaraChannelBuoyData",dir=fdir)
	CSM=getCalcurData("Environ","CapeSanMartinBuoyData",dir=fdir)
	PtSanLuis=getCalcurData("Environ","PtSanLuisDailySST",dir=fdir)
	PtSanLuis$Date=as.Date(as.character(PtSanLuis$Date),format="%m/%d/%Y")
	Nsst=getCalcurData("Environ","NorthernCaliforniaSSTData",dir=fdir)
#  PtArg=sqlFetch(connection,"PtArguelloBuoyData")
#  PtSm=sqlFetch(connection,"PtSantaMariaBuoyData")
#  ESB=sqlFetch(connection,"EastSantaBarbaraChannelBuoyData")
#  WSB=sqlFetch(connection,"WestSantaBarbaraChannelBuoyData")
#  PtSanLuis=sqlFetch(connection,"PtSanLuisDailySST")
#  Nsst=sqlFetch(connection,"NorthernCaliforniaSSTData")
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
  anomalies.bylocation[,,6]=create.anomalies(CSM,yrange,average.years)
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
	 xx=saveCalcurData(PtArguelloSSTAnomalies,"Environ","PtSanLuisSSTAnomalies",dir=fdir)
	 xx=saveCalcurData(PtSanLuisSSTAnomalies,"Environ","PtArguelloSSTAnomalies",dir=fdir)
	 xx=saveCalcurData( PtSantaMariaSSTAnomalies,"Environ","PtSantaMariaSSTAnomalies",dir=fdir)
	 xx=saveCalcurData(WestSantaBarbaraSSTAnomalies,"Environ","WestSantaBarbaraSSTAnomalies",dir=fdir)
	 xx=saveCalcurData(EastSantaBarbaraSSTAnomalies,"Environ","EastSantaBarbaraSSTAnomalies",dir=fdir)
	 xx=saveCalcurData(MontereyBaySSTAnomalies,"Environ","MontereyBaySSTAnomalies",dir=fdir)
	 xx=saveCalcurData(CSMSSTAnomalies,"Environ","CapeSanMartinSSTAnomalies",dir=fdir)
	 xx=saveCalcurData(PtReyesSSTAnomalies,"Environ","PtReyesSSTAnomalies",dir=fdir)
#	 xx=sqlDrop(connection,"PtSanLuisSSTAnomalies",errors=FALSE)
#	 xx=sqlDrop(connection,"PtArguelloSSTAnomalies",errors=FALSE)
#	 xx=sqlDrop(connection,"PtSantaMariaSSTAnomalies",errors=FALSE)
#	 xx=sqlDrop(connection,"WestSantaBarbaraSSTAnomalies",errors=FALSE)
#	 xx=sqlDrop(connection,"EastSantaBarbaraSSTAnomalies",errors=FALSE)
#	 xx=sqlDrop(connection,"MontereyBaySSTAnomalies",errors=FALSE)
#	 xx=sqlDrop(connection,"CapeSanMartinSSTAnomalies",errors=FALSE)
#	 xx=sqlDrop(connection,"PtReyesSSTAnomalies",errors=FALSE)
#	 xx=sqlSave(connection,PtArguelloSSTAnomalies,tablename="PtArguelloSSTAnomalies",append=FALSE,rownames=FALSE)
#	 xx=sqlSave(connection,PtSanLuisSSTAnomalies,tablename="PtSanLuisSSTAnomalies",append=FALSE,rownames=FALSE)
#	 xx=sqlSave(connection,PtSantaMariaSSTAnomalies,tablename="PtSantaMariaSSTAnomalies",append=FALSE,rownames=FALSE)
#	 xx=sqlSave(connection,EastSantaBarbaraSSTAnomalies,tablename="EastSantaBarbaraSSTAnomalies",append=FALSE,rownames=FALSE)
#	 xx=sqlSave(connection,MontereyBaySSTAnomalies,tablename="MontereyBaySSTAnomalies",append=FALSE,rownames=FALSE)
#	 xx=sqlSave(connection,CSMSSTAnomalies,tablename="CapeSanMartinSSTAnomalies",append=FALSE,rownames=FALSE)
#	 xx=sqlSave(connection,PtReyesSSTAnomalies,tablename="PtReyesSSTAnomalies",append=FALSE,rownames=FALSE)
 }
#  odbcClose(connection)
  return(anomalies.bylocation)
}
