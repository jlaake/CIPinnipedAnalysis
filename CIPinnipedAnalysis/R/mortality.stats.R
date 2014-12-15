#' Early pup mortality estimation
#' 
#' Creates a survivorship curve from dead pup count surveys.  The curve extends
#' from birth to either Aug or Sept-Oct depending on the extent of the
#' mortality surveys which has varied across the years.
#' 
#' The script EarlyPupMortality.r runs the functions to construct the
#' ZcEarlPupMortality and CuEarlyPupMortality tables in
#' CIPinnipedCensusQuery.mdb.  This function uses the ZcProduction and
#' CuProduction tables created by PupProduction scrupt and should be
#' run after it is updated with any new data. The CU values use pre-defined correction factors
#' whereas the ZC values use correction factor constructed from tagging data on dead pups. (see popan.cf)
#' 
#' @usage CuMortalityStats(PreLiveCountCf=1.33,PostLiveCountCf=1.25,fdir1,fdir2,years=NULL)
#'        ZcMortalityStats(year,dead,pups)
#' @param PreLiveCountCf Multiplicative correction factor for observed
#'   mortality prior to live count to account for dead pups that were missed
#'   and decomposed or were buried
#' @param PostLiveCountCf Multiplicative correction factor for observed
#'   mortality after live count
#' @param fdir1 database directory for CIPinnipedCensusQuery
#' @param fdir2 database directory for CIPinnipedCensusMaster
#' @param years vector of years that should be selected if not all years
#' @param year year for mortality stats
#' @param dead is a dataframe of corrected dead pup counts across occasions
#' @param pups is the estimated number of pups produced in the mortality count area
#' @export CuMortalityStats ZcMortalityStats
#' @aliases CuMortalityStats ZcMortalityStats
#' @return dataframe containing year,SurveyDate,Pups(number of pups alive at
#'   time),ObservedDead (number observed to die in that interval)
#'   AdjustedDead(corrected count that died),MortalityRate(mortality rate
#'   during the interval between surveys),DailyMortalityRate (mortality rate on
#'   a constant daily rate), CumS (cummulative survival to that time).
#' @author Jeff Laake
CuMortalityStats <-function(PreLiveCountCf=1.33,PostLiveCountCf=1.25,fdir1,fdir2,years=NULL)
{
Castle.Rock.Areas=c("ECR","WCR","CAS")
#
# Read in Cu dead count summary table from the ACCESS file CIPinnipedCensusQuery.mdb
#

  Production=getCalcurData("CIPquery","CuProduction",dir=fdir1)
  CUDates=getCalcurData("CIPCensus","CU Survey Dates",dir=fdir2)
  CUDates$days=as.double(difftime(CUDates$Date,strptime(paste(c(6),c(15),CUDates$Year,sep="/"), "%m/%d/%Y")))
  Production=Production[Production$Area=="Mainland",]
  dead=construct.dead("SMI","Cu",dir=fdir2)
  dead=dead[!dead$MatchingLiveArea%in%Castle.Rock.Areas,]
  dead=dead[dead$DeadPupArea!="N",]
  dead2=with(dead,data.frame(Year=Year,Month=Month,Day=Day,Survey=SurveyNumber,Area=MatchingLiveArea,count=count,type="Dead"))
  dead2=cbind(dead2,Days=as.double(difftime(strptime(paste(dead2$Month,dead2$Day,dead2$Year,sep="/"), "%m/%d/%Y"),strptime(paste(c(6),c(15),dead2$Year,sep="/"), "%m/%d/%Y"))))
  dead2$Days=floor(dead2$Days)
  dead2$Area=factor(dead2$Area)
  dead2=dead2[order(dead2$Area,dead2$Days),]
  counts= tapply(dead2$count, list(as.factor(dead2$Year),factor(dead2$Survey,levels=1:max(CUDates["Survey number"]))), sum)
#
# Cu survey dates and days calculations
#
  days=tapply(CUDates$days,list(CUDates$Year,CUDates[,"Survey number"]),mean) 
  counts=counts[as.numeric(row.names(counts))>1996,]    # restrict analysis to surveys from 1997 and above
  days=days[as.numeric(row.names(days))>1996,]    # restrict analysis to surveys from 1997 and above
  counts[!is.na(days)&is.na(counts)]=0   # assumes that if a survey was done and no cu recorded that none died
  live.dates=as.POSIXct(Production$LiveCountDate[Production$Year>1996])
  mortality.matrix=NULL
  survey.dates=NULL
  for (i in 1:dim(counts)[1])
  {
    year=as.numeric(row.names(counts)[i])
    pups=Production$PupProduction[Production$Year==year]
    pups.produced=pups
    for (j in 1:dim(counts)[2])
    {
       if(!is.na(counts[i,j]))
       {
          survey.date=as.POSIXct(paste(year,"-06-15",sep=""))+floor(days[i,j]+.5)*24*3600
          survey.dates=c(survey.dates,survey.date)
          if(survey.date<=live.dates[i] | is.na(counts[i,2]))
          {
             dead=counts[i,j]*PreLiveCountCf
          }
          else
          {
             if(survey.dates[length(survey.dates)-1]<live.dates[i])
             {
                DeadAtLiveSurvey=Production$DeadInDeadSampleArea[Production$Year==year]
                DeadPostLiveSurvey=sum(counts[i,1:j])-DeadAtLiveSurvey
                DeadPreLiveSurvey=counts[i,j]-DeadPostLiveSurvey
                dead=DeadPostLiveSurvey*PostLiveCountCf + DeadPreLiveSurvey*PreLiveCountCf
             }
             else
                dead=counts[i,j]*PostLiveCountCf
          }
          mr=dead/pups
          cat(year,i,j,days[i,j],mr,"\n")
          if(j==1)
             dmr=1-(1-mr)^(1/days[i,1])
          else
             dmr=1-(1-mr)^(1/(days[i,j]-days[i,j-1]))
          mortality.matrix=rbind(mortality.matrix,c(year,pups,counts[i,j],dead,mr,dmr,cumr=(pups-dead)/pups.produced))
          pups=pups-dead
       }
    }
 }
 df=data.frame(mortality.matrix)
 names(df)=c("Year","Pups","ObservedDead","AdjustedDead","MortalityRate","DailyMortalityRate","CumS")
 # not sure what this is doing since it adds 0 days; must coerce
 df$SurveyDate=survey.dates+as.POSIXct("2006-11-15")-as.numeric(as.POSIXct("2006-11-15"))
return(df)
}
ZcMortalityStats=function(year,dead,pups)
{
#   total across occasions
	cumdead=with(dead,tapply(cumdead,Occasion,sum))
	estdead=with(dead,tapply(estimate.dead,Occasion,sum))
#   compute days from 15 June for mortality rate
	days=with(dead,tapply(daysfrom1July,Occasion,mean))+16
#   smooth estimates if decrease
	if(length(estdead)>=3)
		for(i in 2:(length(estdead)-1))
			if(estdead[i]<estdead[i-1]) estdead[i]=mean(estdead[c(i-1,i+1)])
	if(length(estdead)>1)
	  if(estdead[length(estdead)]<estdead[length(estdead)-1])estdead[length(estdead)-1]=mean(c(estdead[length(estdead)],estdead[length(estdead)-2]))            
#   compute mortality rates
	dead=diff(c(0,estdead))
	mr=dead/pups
	daydiff=diff(c(0,days))
	dmr=1-(1-mr)^(1/(daydiff))
	return(data.frame(Year=rep(year,length(estdead)),Numberdead=dead,DaysFrom15June=days,MortalityRate=mr,DailyMortalityRate=dmr,CumS=1-cumsum(mr)))
}

