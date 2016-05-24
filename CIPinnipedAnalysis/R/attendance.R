#' Create daily attendance encounter history for female sea lions
#' 
#' For a given year and day range (firstday,lastday), finds females that have been seen with a pup in July 
#' in a set of areas (e.g., LTC,LNC,EACS,OCV) and constructs a daily encounter history of seen (1) and not seen(0)
#' to model attendance patterns. 
#' 
#' @param resights dataframe of zalpohus resights from Alive table in database
#' @param year 4 digit year to extract
#' @param firstday numeric value of first day in May to use
#' @param lastday numeric value of last day in July to use
#' @param areas area codes to filter resightings
#' @return list with a vector encounter histories (ch), vector of days when sighting took place (seenDays) and vector of all days (Days). Length of ch will be max of Days.
#' @export
#' @author Jeff Laake
#' @examples 
#' resights=getCalcurData(db="Zc", tbl="Alive")
#' get_attendance(resights,year=2006,firstday=5)
#' 
get_attendance=function(resights,year=2006,firstday=20,lastday=25,areas=c("OCV","EACS","LNC","LTC"))
{
	pupyear=year
	if(year<2000)pupyear=pupyear-1900
	resights=resights[resights$pupyear==pupyear&(as.POSIXlt(resights$sitedate)$mon+1)%in%5:7,]
	resights=resights[as.POSIXlt(resights$sitedate)>=strptime(paste("05/",formatC(firstday,width=2,flag=0),"/",year,sep=""), "%m/%d/%Y"),]
	resights=resights[resights$sitedate<=strptime(paste("07/",formatC(lastday,width=2,flag=0),"/",year,sep=""), "%m/%d/%Y"),]
    resights=resights[resights$code%in%areas,]
	wp=table(resights$brand,resights$withpup)
	pupid=rownames(wp[wp[,"p"]>0 | wp[,"P"]>0 | wp[,"y"]>0 | wp[,"Y"]>0,])
	idpupinJune=resights$brand[(as.POSIXlt(resights$sitedate)$mon+1)==6 & resights$withpup%in%c("p","P","y","Y")]
	idpupinJuly=resights$brand[(as.POSIXlt(resights$sitedate)$mon+1)==7 & resights$withpup%in%c("p","P","y","Y")]
	pupid=pupid[pupid%in%idpupinJuly&pupid%in%idpupinJune] 
	att=droplevels(resights[resights$brand%in%pupid,])
	seenDates=sort(unique(att$sitedate))
	seenDays=as.numeric((seenDates-seenDates[1])/(60*60*24))+1
	dates=strptime(c(paste("05/",formatC(firstday:30,width=2,flag=0),"/",year,sep=""),paste("06/",formatC(1:30,width=2,flag=0),"/",year,sep=""),paste("07/",formatC(1:lastday,width=2,flag=0),"/",year,sep="")), "%m/%d/%Y")
	date=factor(format(att$sitedate),levels=format(dates))
	attch=with(att,table(brand,date))
	attch[attch>0]=1
	attch=attch[rowSums(attch)>2,]
	attch=apply(attch,1,paste,collapse="")
	return(list(ch=attch,seenDays=seenDays,Days=1:length(dates),brand=names(attch)))
}


