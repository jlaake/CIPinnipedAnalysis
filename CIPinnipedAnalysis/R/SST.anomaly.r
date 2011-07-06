

#' Create SST monthly anomaly matrices and dataframes
#' :\code{create.anomalies} computes a matrix of year-month anomalies and
#' \code{create.anomaly.table} converts the matrix into a dataframe of
#' anomalies.
#' 
#' Creates daily means from hourly observations.  Averages daily means within
#' month.  Computes average monthly value from a defined set of years
#' (\code{average.years}) and subtracts the monthly means from the observed
#' year-month matrix to produce anomaly matrix.
#' 
#' @aliases create.anomalies create.anomaly.table
#' @param x for \code{create.anomalies} a dataframe containing WTMP (sea
#'   surface temp) and either a date or 3 separate fields MM,DD,YYYY for
#'   numeric month, day and year respectively; for \code{create.anomaly.table}
#'   x is a matrix that was created by \code{create.anomalies}.
#' @param yrange range of years to construct matrix of anomalies
#' @param average.years vector of years to be included in the monthly mean
#'   calculations the monthly means are subtracted from the specific year-month
#'   values to create the anomaly
#' @export
#' @return for \code{create.anomalies} a matrix of monthly anomalies with years
#'   as rows and months as columns
#' 
#' for \code{create.anomaly.table } a dataframe with fields: Year - 4 digit
#'   numeric, Month - character, midDate - date 15th of the month, and
#'   SSTAnomaly - SST anomaly value.
#' @author Jeff Laake
#' @seealso \code{\link{create.SST.anomalies}}
create.anomalies=function(x,yrange,average.years)
{
#  Function to create a matrix of anomalies.
#  Arguments:
#    x - dataframe containing the following fields:
#          WTMP - Sea surface temp
#          date - Date of measurement
#            or
#                MM    - numeric month
#                DD    - numeric day
#                YYYY  - 4-digit year
#   yrange - range of years for table
#   average.years - years to be used to calculate long-term mean for anomaly measure
#
#  Value: matrix of monthly anomalies with years as rows and months as columns        
#################################################################################
#    Exclude any invalid temps
     x=x[!x$WTMP>=99,]
#    If there is no date field, create one
     if(is.null(x$date))
       x$date=as.Date(paste(formatC(x$MM,width=2,flag=0),"/",formatC(x$DD,width=2,flag=0),"/",x$YYYY,sep=""),format="%m/%d/%Y")
#    Compute daily mean temps
     DailyMeans=tapply(x$WTMP,x$date,mean)
#    Create a dataframe with daily means
     dates=as.Date(names(DailyMeans),format="%Y-%m-%d")
     x=data.frame(Month=as.POSIXlt(dates)$mon,Year=factor(as.POSIXlt(dates)$year+1900,levels=yrange),SST=DailyMeans)
#    Create monthly means from daily means for the set of years defined in average.years
     year.monthly.means=with(x[x$Year%in%average.years,],tapply(SST,list(Month,Year),mean))
#    Compute mean for each month over the set of years
     monthly.means=rowMeans(year.monthly.means,na.rm=TRUE)
#    Compute year/month means for all years and subtract off monthly means to create anomaly matrix and return it
     year.monthly.means=with(x,tapply(SST,list(Month,Year),mean))
     anomalies=t(year.monthly.means-monthly.means)
     return(anomalies)
}
create.anomaly.table=function(x)
{
##################################################################################
# Function to create a dataframe from a matrix of anomalies.  
#   Exclude any months with NaN or NA values.
#
#  Arguments:
#      x - matrix of anomalies with years as rows and months as columns
#
#  Value: dataframe with fields 
#     Year  - 4 digit numeric
#     Month - character
#     midDate - date 15th of the month
#     SSTAnomaly - SST anomaly value
###################################################################################
   yrange=as.numeric(row.names(x))
   dates=paste(rep(yrange,each=12),rep(paste("-",c(rep(0,9),"","",""),1:12,"-15",sep=""),max(yrange)-min(yrange)+1),sep="")
   x=data.frame(Year=rep(yrange,each=12), Month=rep(c("Jan","Feb","Mar","Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"),max(yrange)-min(yrange)+1),
                      midDate=as.character(dates),SSTAnomaly=as.vector(t(x)))
   x=x[!is.nan(x$SSTAnomaly)&!is.na(x$SSTAnomaly),]
   x$midDate=as.character(x$midDate)
   return(x)
}
