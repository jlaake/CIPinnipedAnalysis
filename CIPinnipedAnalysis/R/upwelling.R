#' Upwelling functions
#' 
#' @param uwi dataframe of upwelling (see uwi33,uwi36,uwi39)
#' @param cumuwi cummulative uwi vector
#' @param main title for plot
#' @param colors vector of colors for cummulative upwelling lines
#' @param day julian day value at which cummulative upwelling is extracted by year
#' @param day2 if not null then cummulative upwelling plotted from day to day2
#' @return cum upwelling list from compute.cumuwi, and vector of cummulative 
#' upwelling values by year from compute.cumtodate
#' @usage compute.cumuwi(uwi)
#'        plot_cumuwi(uwi,main="",colors=NULL)
#'        compute.cumtodate(cumuwi,day=NULL)
#'        plot_cumuwi.day(uwi,day=NULL,day2=NULL,main="")
#' @aliases compute.cumuwi plot_cumuwi compute.cumtodate plot_cumuwi.day
#' @author Jeff Laake
#' @export compute.cumuwi plot_cumuwi compute.cumtodate plot_cumuwi.day
#' @examples
#'{
#'  data(uwi33)
#'  data(uwi36)
#'  data(uwi39)
#'	par(mfrow=c(3,1))
#'	colors=rep("grey",2016-1967+1)
#'	colors[c(2012,2014:2015)-1967+1]="green"
#'	plot_cumuwi(uwi39,main="39N",colors=colors)	
#'	plot_cumuwi(uwi36,main="36N",colors=colors)	
#'	plot_cumuwi(uwi33,main="33N",colors=colors)	
#'	
#'	par(mfrow=c(3,1))
#'	plot_cumuwi.day(uwi39,main="39N End of Year")	
#'	plot_cumuwi.day(uwi36,main="36N End of Year")	
#'	plot_cumuwi.day(uwi33,main="33N End of Year")	
#'	
#'	par(mfrow=c(3,1))
#'	plot_cumuwi.day(uwi39,day=151,main="39N to 1 June")	
#'	plot_cumuwi.day(uwi36,day=151,main="36N to 1 June")	
#'	plot_cumuwi.day(uwi33,day=151,main="33N to 1 June")	
#'	
#'	par(mfrow=c(3,1))
#'	plot_cumuwi.day(uwi39,day=90,day2=151,main="39N 1 April to 1 June")	
#'	plot_cumuwi.day(uwi36,day=90,day2=151,main="36N 1 April to 1 June")	
#'	plot_cumuwi.day(uwi33,day=90,day2=151,main="33N 1 April to 1 June")	
#'	
#'}
compute.cumuwi=function(uwi) 
{
	uwi$year=factor(substr(uwi$YYYYMMD,1,4))
	uwi[uwi==-9999]=0
	cumuwi=tapply(uwi$Index,uwi$year,cumsum)
	cumuwi
}
plot_cumuwi=function(uwi,main="",colors=NULL)
{
	cumuwi=compute.cumuwi(uwi)
	if(is.null(colors))colors=rep("black",length(cumuwi))
	plot(1:365,cumuwi[[1]],type="l",ylim=c(min(sapply(cumuwi,min)),max(sapply(cumuwi,max))),xlab="Julian Day",ylab="Cummulative upwelling",main=main,col=colors[1])
	for(i in 2:length(cumuwi))
		lines(1:length(cumuwi[[i]]),cumuwi[[i]],type="l",col=colors[i])
	abline(h=0,col="red",lwd=4)
	invisible()
}	
compute.cumtodate=function(cumuwi,day=NULL)
{
	if(is.null(day)) 
		return(sapply(cumuwi,function(x) x[length(x)]))
	else
		return(sapply(cumuwi,function(x) x[day]))
}	  
plot_cumuwi.day=function(uwi,day=NULL,day2=NULL,main="")
{
	cumuwi=compute.cumuwi(uwi)
	cumtoday=compute.cumtodate(cumuwi,day)
	if(!is.null(day2))
		cumtoday=compute.cumtodate(cumuwi,day2)-cumtoday
	plot(as.numeric(names(cumuwi)),cumtoday,xlab="Year",ylab="Cummulative upwelling",main=main)    
}



