#' Pup production stats
#' 
#' Computes pup production stats for Zalophus or Callorhinus at San Miguel or
#' Castle Rock and creates production tables in the ACCESS database.
#' 
#' Dead pup surveys are conducted typically on 3 occasions during July.  Dead
#' pups are either stacked or tagged so they are not recounted.  A live count
#' is also conducted in late July and it often occurs between two of the dead
#' pup surveys.  This function totals up the dead pups, applies an average
#' multiplicative correction factor to account for pups that died and were
#' never counted, and then interpolates the number that died up to the time of
#' the live count.  The dead pup surveys are only conducted in a portion of the
#' rookery at SMI so the estimated number of dead pups in the dead pup sample
#' area and the live count in that same area are used to create a pre-census
#' mortality rate.  That rate is then used to scale the total live count to
#' compute the number of pups produced (born). The computed values are output
#' in the production table for each species.
#' 
#' \code{\link{do.pup.production}} is the driver function which calls this
#' function for Zc for SMI and Castle Rock and then combines the result in the
#' table ZcProduction which is stored in CIPinnipedCensusQuery.mdb.  Then it
#' does the same thing for Cu and creates CuProduction.
#' 
#' @param island the character name for the island "SMI" or "SNI"
#' @param mainland If TRUE, computation is for San Miguel mainland; otherwise,
#'   Castle Rock
#' @param years vector of years to select or NULL if all to be used
#' @param PreLiveCountCf Multiplicative correction factor for observed
#'   mortality prior to live count to account for dead pups that were missed
#'   and decomposed or were buried
#' @param species either "Zc" for Zalophus or "Cu" for Callorhinus
#' @param dir database directory for CensusMaster
#' @return dataframe containing results
#'   LiveCountDate,Area,Year,LiveInDeadSampleArea,DeadInDeadSampleArea,AdjustedDeadInDeadSampleArea,
#'   MortalityRateAtLiveCount,TotalLiveCountByYear,PupProduction
#' @author Jeff Laake
production.stats <-
function(island="SMI",mainland=TRUE,years=NULL,PreLiveCountCf=1.33,species="Zc",dir=NULL)
{
# Computes pup production stats for San Miguel or Castle Rock
# Define a function to create a dataframe of cummulative counts of dead pups with dates
#  for use in extrapolating number dead at time of live count which typically falls between
#  dead pup counts and thus the need for interpolation.
#  The data frame contains the dates (Days since 15 June) and cummulative count of dead pups up to that date.
summarize.dead=function(x)
{
  counts=sapply(split(x$count,x$SurveyNumber),sum,na.rm=TRUE)
  days=sapply(split(x$Days,x$SurveyNumber),mean,na.rm=TRUE)
  if(length(counts)==0)
     z=data.frame(CumCount=0,Days=days)
  else
  {
     if(length(counts)!=length(days))
     {
       sink()
       stop(paste("\nFor Area = ",x$Area[1]," there are ",length(counts)," counts, and ",length(days)," surveys",sep=""))
     }
     if(any(days!=sort(days)))
     {
       sink()
       stop(paste("\nFor Area = ",x$Area[1]," average of survey dates are not in survey number order",sep=""))     
     }
     z=data.frame(CumCount=cumsum(counts),Days=days)
  }
  return(z)
}
# Start analysis and give title
if(island=="SNI" | mainland)
   cat(paste("\n ", island, "Production Stats for ",species," on mainland\n"))
else
   cat("\nProduction Stats for ",species," on Castle Rock\n")
#
# Read in ZC live and dead count tables from the ACCESS file SMICensusQuery.mdb
#
if(species=="Zc")
{
  live=construct.live(island,"Zc",dir)
  dead=construct.dead(island,"Zc",dir)
}
else
{
  live=construct.live("SMI","Cu",dir)
  dead=construct.dead("SMI","Cu",dir)
  CUDates=getCalcurData("CIPCensus","CU Survey Dates",dir=dir)
#  CUDates=sqlFetch(connection,"CU Survey Dates")
  CUDates$days=as.double(difftime(CUDates$Date,strptime(paste(c(6),c(15),CUDates$Year,sep="/"), "%m/%d/%Y")))
}
if(!is.null(years))
{
   live=live[live$Year%in%years,]
   dead=dead[dead$Year%in%years,]
}   
if(any(is.na(live$Year))) cat("\nWarning: record with missing year in Live file\n")
live=live[!is.na(live$Year),]
if(any(is.na(dead$Year))) cat("\nWarning: record with missing year in Dead file\n")
dead=dead[!is.na(dead$Year),]
#
#  Can be run for either SMI or Castle Rock
#
Castle.Rock.Areas=c("ECR","WCR","CAS")
if(island=="SMI")
{
  if(mainland)
  {
     live=live[!live$Area%in%Castle.Rock.Areas,]
     dead=dead[!dead$MatchingLiveArea%in%Castle.Rock.Areas,]
  }
  else
  {
     live=live[live$Area%in%Castle.Rock.Areas,]
     dead=dead[dead$MatchingLiveArea%in%Castle.Rock.Areas,]
  }
}
#
# Stop if missing a DeadPupArea value
#
if(any(is.na(live$DeadPupArea)))
{
   cat("\nMissing DeadPupArea value in live file:\n")
   print(live[is.na(live$DeadPupArea),] )
   sink()
   stop("\nProduction table not created\n")
}
if(any(is.na(dead$DeadPupArea)))
{
   cat("\nMissing DeadPupArea value in dead file:\n")
   print(dead[is.na(dead$DeadPupArea),] )
   sink()
   stop("\nProduction table not created\n")
}
#
# Exclude areas that were not included in the dead pup sampling
#
live2=with(live[live$DeadPupArea!="N",],data.frame(Year=Year,Month=Month,Day=Day,Area=Area,count=AvgCount,type="Live"))
dead2=with(dead[dead$DeadPupArea!="N",],data.frame(Year=Year,Month=Month,Day=Day,Area=MatchingLiveArea,SurveyNumber=SurveyNumber,count=count,type="Dead"))
# Exclude months before June
dead2=dead2[dead2$Month>= 6,]
#
# Check to make sure that each count was assigned a value for DeadPupArea and MatchingLiveArea in the DeadPupSamplesArea Access table
#
if(any(levels(dead2$Area)==""))
{
  cat("\nFollowing records are missing the MatchingLiveArea and DeadPupArea field\n")
  dead2[dead2$Area=="",]
  sink()
  stop("\n",species," Production table not constructed")
}
#
#  Recreate the Area factor to drop levels that aren't used
#
live2$Area=factor(live2$Area)
#
#  Add a field to live and dead count which is the number of days passed 15 June of each year
#
live=cbind(live,Days=as.double(difftime(strptime(paste(live$Month,live$Day,live$Year,sep="/"), "%m/%d/%Y"),strptime(paste(c(6),c(15),live$Year,sep="/"), "%m/%d/%Y"))))
live2=cbind(live2,Days=as.double(difftime(strptime(paste(live2$Month,live2$Day,live2$Year,sep="/"), "%m/%d/%Y"),strptime(paste(c(6),c(15),live2$Year,sep="/"), "%m/%d/%Y"))))
dead2=cbind(dead2,Days=as.double(difftime(strptime(paste(dead2$Month,dead2$Day,dead2$Year,sep="/"), "%m/%d/%Y"),strptime(paste(c(6),c(15),dead2$Year,sep="/"), "%m/%d/%Y"))))
live2$Days=floor(live2$Days)
dead2$Days=floor(dead2$Days)
#
#  Apportion the dead pup count from 2002 that was combined for NWP/NEP/NWC. This is only 77 pups and it is
#  apportioned based on amount of mortality in that area from the other surveys that year
#
if(island=="SMI" & mainland &species=="Zc")
{
   splitcounts=tapply(dead2$count[dead2$Year==2002&dead2$Area%in%c("NWP","NEP","NWC")],dead2$Area[dead2$Year==2002&dead2$Area%in%c("NWP","NEP","NWC")],mean)
   splitcounts=splitcounts[!is.na(splitcounts)]
   splitcounts=splitcounts/sum(splitcounts)
   splitcounts=splitcounts*dead2$count[dead2$Area=="NEP/NWC/NWP"&dead2$Year==2002]
   dead2=rbind(dead2,dead2[dead2$Area=="NEP/NWC/NWP"&dead2$Year==2002,],dead2[dead2$Area=="NEP/NWC/NWP"&dead2$Year==2002,])
   dead2$count[dead2$Area=="NEP/NWC/NWP"&dead2$Year==2002]=splitcounts
   dead2$Area[dead2$Area=="NEP/NWC/NWP"&dead2$Year==2002]=names(splitcounts)
}
dead2$Area=factor(dead2$Area)
dead2=dead2[order(dead2$Area,dead2$Days),]
#
# Create a data frame that has a record for each unique area and year for the live count.
# In some instances the live count occured over more than one day, so the counts are
# summed and the days from 15 June are averaged.
#
LiveByYearandArea=expand.grid(Year=unique(live2$Year),Area=unique(live2$Area))
LiveByYearandArea=LiveByYearandArea[order(LiveByYearandArea$Area,LiveByYearandArea$Year),]
LiveByYearandArea$Count= sapply(split(live2$count,list(live2$Year,live2$Area)),sum)
LiveByYearandArea$Days= sapply(split(live2$Days,list(live2$Year,live2$Area)),mean)
LiveByYearandArea=LiveByYearandArea[!is.nan(LiveByYearandArea$Days),]
#
# Using summarize.dead function create a similar structure for dead pups except that it is a list
# with one entry per year-area.  Each list element contains a data frame with
# the dates (Days since 15 June) and cummulative count of dead pups up to that date.
#
DeadByYearandArea=lapply(split(subset(dead2,select=c("count","Days","SurveyNumber","Area")),list(dead2$Year,dead2$Area),drop=TRUE),summarize.dead)
#
# It is essential that the live and dead dataframes match up in length and have the same
# structure for year-area.  If the dimensions are not matching say so.
#
if(dim(LiveByYearandArea)[1]!=length(DeadByYearandArea))
{
  cat("\nThe length of the live counts does not match length of dead counts\n")
  cat("\nLength of live = ",dim(LiveByYearandArea)[1], "Length of dead = ",length(DeadByYearandArea),"\n")
}
#
#  Next check to make sure that for every live count year-area there is a dead count
#  and vice-versa.
#
livenames=paste(live2$Year,live2$Area,sep=".")
missing.dead.areas=livenames[!livenames%in%names(DeadByYearandArea)]
if(length(missing.dead.areas)!=0)
{
  cat("The following live pup counts are missing a matching pup mortality count in this year and area:",paste("\n",missing.dead.areas,sep=""),"\n")
  cat("Mortality for these areas is assumed to be 0.\n")
}
missing.live.areas=names(DeadByYearandArea)[!names(DeadByYearandArea)%in%livenames]
if(length(missing.live.areas)!=0)
{
  cat("The following dead pup counts are missing a matching live pup count in this year and area:", paste("\n",missing.live.areas,sep=""),"\n")
  sink()
  stop("\n",species," Production table not constructed")
}
#
# If no errors, then for each live area interpolate cumulative mortality at date of live count
# and store as Dead.at.live.count.  If all of the dead pup counts are after the live pup count
# this is an error and should not happen.
#
livenames=paste(LiveByYearandArea$Year,LiveByYearandArea$Area,sep=".")
LiveByYearandArea$Dead.at.live.count=0
for (i in 1:length(livenames))
{        
     if(is.null(DeadByYearandArea[[livenames[i]]]$Days))
        dead.at.livecount=0
     else
     {
           xyear=as.numeric(substr(livenames[i],1,4))
           if(species=="Zc" | xyear < 1997 | island=="SNI")
           {
              if((island=="SNI" | species=="Cu") & length(DeadByYearandArea[[livenames[i]]]$Days)==1)
                 dead.at.livecount=DeadByYearandArea[[livenames[i]]]$CumCount
              else
              {
                x=c(0,DeadByYearandArea[[livenames[i]]]$Days)
                y=c(0,DeadByYearandArea[[livenames[i]]]$CumCount)
                if(length(x)!=length(y))
                  cat(paste("\nWarning: For ",livenames[i]," length of average survey dates and counts do not match"))
                if(length(x)!=length(unique(x)))
                  cat(paste("\nWarning: For ",livenames[i]," some average survey dates are not unique"))
			    if(LiveByYearandArea$Days[i]>x[length(x)])
			    {
				    cat(paste("For",livenames[i],"shifting time of livecount from day ",LiveByYearandArea$Days[i],"to ",x[length(x)]))
				    LiveByYearandArea$Days[i]=x[length(x)]
			    }
			    dead.at.livecount=approx(x=x,y=y,xout=LiveByYearandArea$Days[i])$y
              }
           }
           else
           {
              if(length(DeadByYearandArea[[livenames[i]]]$Days)==1)
                dead.at.livecount=DeadByYearandArea[[livenames[i]]]$CumCount
              else
              {
                 SurveyNum=as.numeric(row.names(DeadByYearandArea[[livenames[[i]]]]))
                 xdays=c(0,sort(CUDates$days[CUDates$Year==xyear&CUDates[,"Survey number"]%in%SurveyNum]))
                 y=c(0,DeadByYearandArea[[livenames[i]]]$CumCount)
                 if(length(xdays)!=length(y))
                   cat(paste("\nWarning: For ",livenames[i]," length of average survey dates and counts do not match"))
                 if(length(xdays)!=length(unique(xdays)))
                   cat(paste("\nWarning: For ",livenames[i]," some average survey dates are not unique"))
			     if(LiveByYearandArea$Days[i]>xdays[length(xdays)])
				 {
					 cat(paste("For",livenames[i],"shifting time of livecount from day ",LiveByYearandArea$Days[i],"to ",xdays[length(xdays)]))
					 LiveByYearandArea$Days[i]=xdays[length(xdays)]
				 }
                 dead.at.livecount=approx(x=xdays,y=y,xout=LiveByYearandArea$Days[i])$y
              }
           } 
     }
	 if(is.na(dead.at.livecount))
           cat(paste("\nProblem with computation of dead at livecount for ",livenames[i]))
     LiveByYearandArea$Dead.at.live.count[i]=dead.at.livecount
}
TotalLiveCountByYear=round(sapply(split(live$AvgCount,live$Year),sum))
Years=as.numeric(names(TotalLiveCountByYear))
mc=getCalcurData("CIPCensus","MissingCounts",dir=dir)
#mc=sqlFetch(connection,"MissingCounts")
mc$Year=factor(mc$Year,levels=Years)
mc=mc[mc$Island==island & mc$Species==species,]
if(island=="SMI")
if(mainland)
	mc=mc[!mc$Area%in%c("ECR","WCR"),]
else
	mc=mc[mc$Area%in%c("ECR","WCR"),]
if(nrow(mc)>0)
{
	mc=sapply(split(mc$Count,mc$Year),sum)
	TotalLiveCountByYear=TotalLiveCountByYear+mc
}	
LiveInDeadSampleArea=round(sapply(split(LiveByYearandArea$Count,LiveByYearandArea$Year),sum))
DeadInDeadSampleArea=round(sapply(split(LiveByYearandArea$Dead.at.live.count,LiveByYearandArea$Year),sum))
AdjustedDeadInDeadSampleArea=round(PreLiveCountCf*DeadInDeadSampleArea)
MortalityRateAtLiveCount=AdjustedDeadInDeadSampleArea/(AdjustedDeadInDeadSampleArea+LiveInDeadSampleArea)
if(mainland | island=="SNI")
   Area = "Mainland"
else
   Area = "CastleRock"
dfl=data.frame(Year=as.numeric(names(TotalLiveCountByYear)),Live=TotalLiveCountByYear)
dfm=data.frame(Year=as.numeric(names(MortalityRateAtLiveCount)),LD=LiveInDeadSampleArea,ADD=AdjustedDeadInDeadSampleArea,DD=DeadInDeadSampleArea,mr=MortalityRateAtLiveCount)
dfp=merge(dfl,dfm,all.x=TRUE)
dfp$mr[is.na(dfp$mr)]=0
dateoflive=floor(tapply(live$Days,live$Year,mean)+.5)*24*3600+as.POSIXct(paste(unique(live$Year),"-06-15",sep=""))
dat=data.frame(Island=island,LiveCountDate=as.character(dateoflive),Area=rep(Area,length(Years)),Year=Years,LiveInDeadSampleArea=round(dfp$LD),DeadInDeadSampleArea=round(dfp$DD),AdjustedDeadInDeadSampleArea=round(dfp$ADD),
            MortalityRateAtLiveCount=round(dfp$mr*100,2),TotalLiveCountByYear=TotalLiveCountByYear,
            PupProduction=round(dfp$Live/(1-dfp$mr)))
}

