#' Constructs dataframes of live or dead pup counts
#' 
#' Extracts data from CIPinnipedMaster database and constructs dataframe of
#' dead pup counts or live pup counts for specified island and species.
#' 
#' \code{construct.live} extracts data from table \code{Zc Cu live pup census}
#' and sums counts for an observer within an area for a year and then averages
#' those counts across observers that counted the same area that year.
#' 
#' \code{construct.dead} extracts data from table \code{Zc Cu dead pup census}
#' and table \code{Zc dead tag initial} and totals dead pups by survey
#' number/date.  It will produce errors if there are more than one survey date
#' for each survey number within an area in a year.
#' 
#' @aliases construct.live construct.dead
#' @param island Either "SMI" or "SNI"
#' @param species Either "Zc" or "Cu"
#' @param dir database directory
#' @return for \code{construct.live} a dataframe containing YearArea, Survey
#'   date, Area, AvgCount, Year, Month, Day, DeadPupArea (Y/N if included in
#'   area sampled for dead pups).
#' 
#' for \code{construct.dead} a dataframe containing YearArea, SurveyNumber,
#'   SurveyDate, Area, count, Year, Month, Day, DeadPupArea (Y/N if included in
#'   area sampled for dead pups), and MatchingLiveArea.
#' @author Jeff Laake
construct.live=function(island,species,dir=NULL)
{
# Function to create the dataframe of live counts; this used to be done with
# query from Access Zc SMI Live Count by area, year. Instead it now uses table Zc Cu live pup census
# and constructs the total counts by averaging over observer counts in the same area-year.  It can
# also be used for SNI counts.
#
	livepups=getCalcurData("CIPCensus","Zc Cu live pup census",dir=dir)
#	livepups=sqlFetch(connection,"Zc Cu live pup census")
	livepups=livepups[!is.na(livepups$Count),]
	livepups$Observer=as.character(livepups$Observer)
	livepups$Observer[is.na(livepups$Observer)]="Unknown"
	livepups$Observer=factor(livepups$Observer)
	livepups[is.na((livepups[,"Survey date"])),][,"Survey date"]=as.Date(paste("1Aug",livepups$Year[is.na(livepups[,"Survey date"])],sep=""),"%d%b%Y")
	livepups=livepups[livepups$Island==island & livepups$Species==species,]
	count.by.observer=tapply(livepups$Count,paste(livepups$Observer,livepups[,"Survey date"],livepups[,"Area code"]),sum,na.rm=TRUE)
	count.df=as.data.frame(do.call("rbind",strsplit(names(count.by.observer)," ")))
	if(ncol(count.df)!=3)
		count.df=count.df[,-3]
	names(count.df)=c("Observer","Survey date","Area code")
	count.df$Count=count.by.observer
	total.count.by.area=tapply(count.df$Count,paste(count.df[,"Survey date"],count.df[,"Area code"]),mean,na.rm=TRUE)
	count.df=as.data.frame(do.call("rbind",strsplit(names(total.count.by.area)," ")))
	names(count.df)=c("Survey date","Area")
	count.df$AvgCount=total.count.by.area
	count.df[,"Survey date"]=as.Date(count.df[,"Survey date"])
	count.df$Year=as.POSIXlt(count.df[,"Survey date"])[[6]]+1900
	count.df$Month=as.POSIXlt(count.df[,"Survey date"])[[5]]+1
	count.df$Day=as.POSIXlt(count.df[,"Survey date"])[[4]]
	count.df$YearArea=paste(count.df$Year,count.df$Area,sep="")
	dps=getCalcurData("CIPCensus","DeadPupSampleAreas",dir=dir)
#	dps=sqlFetch(connection,"DeadPupSampleAreas")
	dps=dps[dps$YearAreaSpecies==species,]
	tt=table(paste(dps$Year,dps$Area))
	if(any(tt>1))
		stop(paste("Following dead pup sample areas are duplicated for species",species,": "),paste(names(tt)[tt>1],sep=", "))
	count.df=merge(count.df,subset(dps,select=c("YearArea","Dead pup sample area")),by="YearArea")
	names(count.df)[length(names(count.df))]="DeadPupArea"
	return(count.df)
}

construct.dead=function(island,species,dir=NULL)
{
# Function to create the dataframe of dead counts; this used to be done with
# query from Access "Zc SMI Total Dead by Area and Survey Number" or
# "Cu SMI Total Dead by Area and Survey Number". Instead it now uses table Zc Cu dead pup census table
# and Zc dead tag initial table and constructs the total counts by adding values
# by survey date-area across both tables. It can also be used for SNI counts.
#
#   Get dead pup sample area table
    dps=getCalcurData("CIPCensus","DeadPupSampleAreas",dir=dir)
#	dps=sqlFetch(connection,"DeadPupSampleAreas")
	dps=dps[dps$YearAreaSpecies==species,]
    tt=table(paste(dps$Year,dps$Area))
    if(any(tt>1))
     stop(paste("Following dead pup sample areas are duplicated for species",species,": "),paste(names(tt)[tt>1],sep=", "))
#   Get dead stacked pups and filter to island and species and only use fullterm pups
#   and exclude any NA survey number
    dead.stack=getCalcurData("CIPCensus","Zc Cu dead pup census",dir=dir)
#	dead.stack=sqlFetch(connection,"Zc Cu dead pup census")
	dead.stack=dead.stack[dead.stack$Island==island&dead.stack$Species==species&
                                  dead.stack$Development%in%c("fullterm","Fullterm")
                                  &!is.na(dead.stack[,"Survey number"])&dead.stack$Disposition!="Tagged",]
    dead.stack[is.na(dead.stack[,"Survey date"]),"Survey date"]=as.Date(paste("08-01-",dead.stack$Year[is.na(dead.stack[,"Survey date"])],sep=""),"%m-%d-%Y")
    dead.stack[,"Area code"]=factor(dead.stack[,"Area code"])
#   Make sure that there is only one survey number per survey date in each area
    xx=(tapply(dead.stack[,"Survey number"],list(dead.stack[,"Survey number"],dead.stack[,"Survey date"],dead.stack[,"Area code"]),length))
    xx=apply(xx,c(2,3),function(x) length(x[!is.na(x)]))
    if(any(xx>1))
    {
      indices=which(xx>1,arr.ind=TRUE)
      stop(paste("In dead stacked, for species",species, "the following survey dates/areas have more than one survey number: "),paste(rownames(xx)[indices[,"row"]]
          ,colnames(xx)[indices[,"col"]],sep="-",collapse=","))
    }
#   Create dataframe of counts
    xx=sapply(split(dead.stack[,"Number dead"],list(dead.stack[,"Survey number"],dead.stack[,"Survey date"],dead.stack[,"Area code"])),sum,na.rm=TRUE)
    xx=xx[xx>0]
    xmat=as.data.frame(do.call("rbind",strsplit(names(xx),"\\.")))
    names(xmat)=c("SurveyNumber","SurveyDate","Area")
    xmat$SurveyNumber=as.numeric(as.character(xmat$SurveyNumber))
    xmat$SurveyDate=as.Date(as.character(xmat$SurveyDate),format="%Y-%m-%d")
    xmat$count=xx
    xmat$Year=as.POSIXlt(xmat$SurveyDate)[[6]]+1900
    if(species=="Zc")
    {
  	   dead.tag=getCalcurData("CIPCensus","Zc dead tag initial",dir=dir)
#	   dead.tag=sqlFetch(connection,"Zc dead tag initial")
       dead.tag=dead.tag[dead.tag$Island==island&(dead.tag$Development=="fullterm" |dead.tag$Development=="Fullterm"),]
       dead.tag[,"Area code"]=factor(dead.tag[,"Area code"])

#   Make sure that there is only one survey number per survey date in each area
       xx=(tapply(dead.tag[,"Survey number"],list(dead.tag[,"Survey number"],dead.tag[,"Survey date"],dead.tag[,"Area code"]),length))
       xx=apply(xx,c(2,3),function(x) length(x[!is.na(x)]))
       if(any(xx>1))
       {
         indices=which(xx>1,arr.ind=TRUE)
         stop("In dead tag, for Zc, the following survey dates/areas have more than one survey number: ",paste(rownames(xx)[indices[,"row"]]
            ,colnames(xx)[indices[,"col"]],sep="-",collapse=","))
       }
#   Create dataframe of counts
       xx=sapply(split(dead.tag[,"Survey number"],list(dead.tag[,"Survey number"],dead.tag[,"Survey date"],dead.tag[,"Area code"])),length)
       xx=xx[xx>0]
       xmat1=as.data.frame(do.call("rbind",strsplit(names(xx),"\\.")))
       names(xmat1)=c("SurveyNumber","SurveyDate","Area")
       xmat1$SurveyNumber=as.numeric(as.character(xmat1$SurveyNumber))
       xmat1$SurveyDate=as.Date(as.character(xmat1$SurveyDate),format="%Y-%m-%d")
       xmat1$count=xx
       xmat1$Year=as.POSIXlt(xmat1$SurveyDate)[[6]]+1900
       xmat=rbind(xmat,xmat1)
    }
    xmat$YearArea=paste(xmat$Year,xmat$Area,sep="")
    xmat$Month=as.POSIXlt(xmat$SurveyDate)[[5]]+1
    xmat$Day=as.POSIXlt(xmat$SurveyDate)[[4]]
    xmat=merge(xmat,subset(dps,select=c("YearArea","Dead pup sample area","MatchingLiveArea")),by="YearArea",all.x=TRUE)
    names(xmat)[length(names(xmat))-1]="DeadPupArea"
    return(xmat)
}
