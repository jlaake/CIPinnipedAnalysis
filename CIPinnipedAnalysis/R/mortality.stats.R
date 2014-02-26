#' Performs numerical computations for early pup mortality estimation
#' 
#' Creates a survivorship curve from dead pup count surveys.  The curve extends
#' from birth to either Aug or Sept-Oct depending on the extent of the
#' mortality surveys which has varied across the years.
#' 
#' \code{\link{do.early.pup.mortality}} is the driver function to construct the
#' ZcEarlPupMortality and CuEarlyPupMortality tables in
#' CIPinnipedCensusQuery.mdb.  This function uses the ZcProduction and
#' CuProduction tables created by \code{\link{do.pup.production}} and should be
#' run after it is updated with any new data.
#' 
#' @param PreLiveCountCf Multiplicative correction factor for observed
#'   mortality prior to live count to account for dead pups that were missed
#'   and decomposed or were buried
#' @param PostLiveCountCf Multiplicative correction factor for observed
#'   mortality after live count
#' @param species either "Zc" for Zalophus or "Cu" for Callorhinus
#' @param island the character name for the island "SMI" or "SNI"
#' @param fdir1 database directory for CIPinnipedCensusQuery
#' @param fdir2 database directory for CIPinnipedCensusMaster
#' @param years vector of years that should be selected if not all years
#' @return dataframe containing year,SurveyDate,Pups(number of pups alive at
#'   time),ObservedDead (number observed to die in that interval)
#'   AdjustedDead(corrected count that died),MortalityRate(mortality rate
#'   during the interval between surveys),DailyMortalityRate (mortality rate on
#'   a constant daily rate), CumS (cummulative survival to that time).
#' @author Jeff Laake
mortality.stats <-function(PreLiveCountCf=1.33,PostLiveCountCf=1.25,species="Zc",island="SMI",fdir1,fdir2,years=NULL)
{
# Read in ZC dead count summary table from the ACCESS file CIPinnipedCensusQuery.mdb
#
Castle.Rock.Areas=c("ECR","WCR","CAS")
if(species=="Zc")
{
  Production=getCalcurData("CIPquery","ZcProduction",dir=fdir1)
# Production=sqlFetch(connection1,"ZcProduction")
  Production=Production[Production$Area=="Mainland"&Production$Island==island,]
  dead=construct.dead(island,species,dir=fdir2)
  dead=dead[dead$DeadPupArea!="N",]
  if(!is.null(years)) dead=dead[dead$Year%in%years,]
  if(island=="SMI") dead=dead[!dead$MatchingLiveArea%in%Castle.Rock.Areas,]
  dead2=with(dead,data.frame(Year=Year,Month=Month,Day=Day,Survey=SurveyNumber,Area=MatchingLiveArea,count=count,type="Dead"))
}
else
#
# Read in Cu dead count summary table from the ACCESS file CIPinnipedCensusQuery.mdb
#
{
  Production=getCalcurData("CIPquery","CuProduction",dir=fdir1)
  CUDates=getCalcurData("CIPCensus","CU Survey Dates",dir=fdir2)
# Production=sqlFetch(connection1,"CuProduction")
# CUDates=sqlFetch(connection2,"CU Survey Dates")
  CUDates$days=as.double(difftime(CUDates$Date,strptime(paste(c(6),c(15),CUDates$Year,sep="/"), "%m/%d/%Y")))
  Production=Production[Production$Area=="Mainland",]
  dead=construct.dead(island,species,dir=fdir2)
  dead=dead[!dead$MatchingLiveArea%in%Castle.Rock.Areas,]
  dead=dead[dead$DeadPupArea!="N",]
  dead2=with(dead,data.frame(Year=Year,Month=Month,Day=Day,Survey=SurveyNumber,Area=MatchingLiveArea,count=count,type="Dead"))
}
#
# Adjust Zc counts for 2002 when areas were combined for some of the surveys.
#
if(species=="Zc"&island=="SMI")
{
   splitcounts=tapply(dead2$count[dead2$Year==2002&dead2$Area%in%c("NWP","NEP","NWC")],dead2$Area[dead2$Year==2002&dead2$Area%in%c("NWP","NEP","NWC")],mean)
   splitcounts=splitcounts[!is.na(splitcounts)]
   splitcounts=splitcounts/sum(splitcounts)
   splitcounts=splitcounts*dead2$count[dead2$Area=="NEP/NWC/NWP"&dead2$Year==2002]
   dead2=rbind(dead2,dead2[dead2$Area=="NEP/NWC/NWP"&dead2$Year==2002,],dead2[dead2$Area=="NEP/NWC/NWP"&dead2$Year==2002,])
   dead2$count[dead2$Area=="NEP/NWC/NWP"&dead2$Year==2002]=splitcounts
   dead2$Area[dead2$Area=="NEP/NWC/NWP"&dead2$Year==2002]=names(splitcounts)
}
dead2=cbind(dead2,Days=as.double(difftime(strptime(paste(dead2$Month,dead2$Day,dead2$Year,sep="/"), "%m/%d/%Y"),strptime(paste(c(6),c(15),dead2$Year,sep="/"), "%m/%d/%Y"))))
dead2$Days=floor(dead2$Days)
dead2$Area=factor(dead2$Area)
dead2=dead2[order(dead2$Area,dead2$Days),]
if(species=="Zc")
   counts=tapply(dead2$count,list(as.factor(dead2$Year),as.factor(dead2$Survey)),sum)
else
   counts= tapply(dead2$count, list(as.factor(dead2$Year),factor(dead2$Survey,levels=1:max(CUDates["Survey number"]))), sum)

#
# Zc survey dates and days calculation
#
if(species=="Zc")
{
  live.dates=as.POSIXct(Production$LiveCountDate[Production$Year>1990])
  days=tapply(dead2$Days,list(as.factor(dead2$Year),as.factor(dead2$Survey)),mean)
  counts=counts[as.numeric(row.names(counts))>1990,]
  days=days[as.numeric(row.names(days))>1990,]
}
else
#
# Cu survey dates and days calculations
#
{
  days=tapply(CUDates$days,list(CUDates$Year,CUDates[,"Survey number"]),mean) 
  counts=counts[as.numeric(row.names(counts))>1996,]    # restrict analysis to surveys from 1997 and above
  days=days[as.numeric(row.names(days))>1996,]    # restrict analysis to surveys from 1997 and above
  counts[!is.na(days)&is.na(counts)]=0   # assumes that if a survey was done and no cu recorded that none died
  live.dates=as.POSIXct(Production$LiveCountDate[Production$Year>1996])
}
mortality.matrix=NULL
survey.dates=NULL
for (i in 1:dim(counts)[1])
{
    year=as.numeric(row.names(counts)[i])
    if(species=="Zc")
       pups=Production$AdjustedDeadInDeadSampleArea[Production$Year==year]+Production$LiveInDeadSampleArea[Production$Year==year]
    else
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
          {
             dmr=1-(1-mr)^(1/days[i,1])
          }
          else
             dmr=1-(1-mr)^(1/(days[i,j]-days[i,j-1]))
          mortality.matrix=rbind(mortality.matrix,c(year,pups,counts[i,j],dead,mr,dmr,cumr=(pups-dead)/pups.produced))
          pups=pups-dead
       }
    }
}
df=data.frame(mortality.matrix)
names(df)=c("Year","Pups","ObservedDead","AdjustedDead","MortalityRate","DailyMortalityRate","CumS")
df$SurveyDate=survey.dates+as.POSIXct("2006-11-15")-as.numeric(as.POSIXct("2006-11-15"))
return(df)
}

