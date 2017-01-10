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
#' run after it is updated with any new data. The ZC and CU values use correction factors 
#' constructed from tagging data on dead pups. (see popan.cf and cu_popan.cf)
#' 
#' @param year year for mortality stats
#' @param dead is a dataframe of corrected dead pup counts across occasions
#' @param pups is the estimated number of pups produced in the mortality count area
#' @export 
#' @return dataframe containing year,Pups(number of pups alive at
#'   time),Numbersead (number estimated to die in that interval)
#'   DaysFrom15June,MortalityRate(mortality rate
#'   during the interval between surveys),DailyMortalityRate (mortality rate on
#'   a constant daily rate), CumS (cummulative survival to that time).
#' @author Jeff Laake
MortalityStats=function(year,dead,pups)
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

