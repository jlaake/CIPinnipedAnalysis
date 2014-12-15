#' Correction factors for dead pup counts
#' 
#' Correction factors for a set of dead pup counts are created for each survey date or a single date from 
#' estimates created from POPAN model (cfdata) on tagging data.
#' 
#' For a given island and year, it calls getdead_ch  and process_ch to create the capture history for the dead pup data and then 
#' uses process.data from RMark to process the capture history data for the POPAN model. If cfyear > 1997 it uses Position and Substrate as group factors and if 
#' cfyear < 1997 it uses Area code which is split into 4 areas: PTS(NEP,NWP,PBP),SCV (SCV,PBS,EAC and WAC), WCV and NWC. The argument cfdata
#' contains the estimates from the POPAN model to create the correction factor. It is created by the function popan.cf which runs a set of POPAN models
#' with the tagging data from year yyyy (1994,1995,1998,2002 for SMI and 2006 for SNI) and creates smiyyyy.popan.results or sniyyyypopan.results with the model results and model averaged values for cfdata. If ndays is null it creates a correction factor for each 
#' survey occasion. If ndays, number of days from 1 July, is not null then it computes with interpolation the number observed dead at that date and
#' its correction factor.  
#' 
#' The function average_cf computes an average correction factor for a set of years with tagging data. This is only useful at present with SMI becuse SNI has only 
#' one year of tagging data but it can be used with SNI as well.  It can use an area (SCV,PTS,WCV,NWC) or lev value (AC,AN,BC,BN) to compute specific correction 
#' factors. Based on the island and year, it knows which set of years to use for the correction factor.
#' 
#' Note: In 1994 the number of untagged and stacked was not recorded on the last occasion; thus the correction factor for the last occasion is not useful.  
#' 
#' @usage correct_dead(island,year)
#'        average_cf(island,year,ndays1=Inf,ndays2,area=NULL,lev=NULL)
#'        compute_cf(BiGross,ndays1=Inf,ndays2)
#'        process_ch(ch,freq=NULL,all=FALSE)
#' @aliases correct_dead average_cf compute_cf process_ch
#' @export correct_dead average_cf compute_cf process_ch
#' @param island ("SMI" or "SNI")
#' @param year four digit numeric year
#' @param BiGross is the estimates of gross immigration into the population of dead pups; it is a list with results from fitted models
#' @param ndays1 number of days from 1 July for beginning of interval; -Inf default which starts at beginning
#' @param ndays2 number of days from 1 July at end of interval (time of pup count)
#' @param area for computing correction factor with tagging data prior to 1998
#' @param lev position on beach and substrate (AC,AN,BC,BN) for computing correction factor with data in 1998 and later
#' @param ch vector of capture histories
#' @param freq frequency of capture history
#' @param all logical indicating whether all quantities should be computed from capture history
#' @return correct_dead returns a list containing a dataframe with corrected dead pup counts by occasion and strata (eg group - area or substrate/position) and a dataframe with a total across strata by occasion. 
#' average_cf returns an average correction factor for a particular range of dates by interpolation and then averaging over years (cfyears) for a particular area code or postion-substrate strata.  
#' @author Jeff Laake 
#' @seealso popan.cf getdead_ch
correct_dead=function(island,year)
{
#   get dead pups for the island and year as capture histories
	x.popan=suppressMessages(getdead_ch(island,year,merge=FALSE))
#   only use first 3 characters of the area code
	x.popan$df[,"Area code"]=factor(substr(as.character(x.popan$df[,"Area code"]),1,3))
#   ignore any 0 records
	x.popan$df=x.popan$df[!x.popan$df$freq==0,]
#   get first capture occasion from the capture history
	chlist=CIPinnipedAnalysis::process_ch(x.popan$df$ch,x.popan$df$freq)
	times=x.popan$days
	time.intervals=as.numeric(diff(times))
#   use RMark code to process the capture history 
	if(year>1997)
	{
		x.proc=process.data(x.popan$df,model="POPAN",time.intervals=time.intervals,groups=c("Area code","Position","Substrate"),begin.time=0)
	}else
		x.proc=process.data(x.popan$df,model="POPAN",time.intervals=time.intervals,groups=c("Area code"),begin.time=0)
#   tally the number dead at first occasion by group	
	numdead= tapply(abs(x.proc$data$freq),list(factor(chlist$first,1:nrow(times)),x.proc$data$group),sum)
	numdead[is.na(numdead)]=0
#   compute cumulative dead across occasions within each group (strata)
	cumdead=as.vector(apply(t(numdead), 1, cumsum))
#   create a dataframe (df) with occasion and strata (either area or area-substrate-position)
	k=length(times)
	gindex=rep(1:nrow(x.proc$group.covariates),each=k)
	df=x.proc$group.covariates[gindex,,drop=FALSE]
	df$occasion=rep(1:k,nrow(x.proc$group.covariates))
	df$key=paste(df[,"Area code"],df$occasion)
    x.popan$daysfrom1July$key=paste(x.popan$daysfrom1July$Area,x.popan$daysfrom1July$Occasion)
	df=merge(df, x.popan$daysfrom1July,by="key")
	df$occasion=NULL
	df[,"Area code"]=NULL
	df$key=NULL
	if(year>1997)
		df=df[order(df$Substrate,df$Position,df$Area,df$Occasion),]
	else
		df=df[order(df$Area,df$Occasion),]
	df$cumdead=cumdead
#
#   compute and apply correction factor to observed data
#
 	df$cf=rep(NA,nrow(df))
	for(i in 1:nrow(df))
	{
		if(year<=1997)
		   df$cf[i]=average_cf(island=island,year=year,ndays2=df$daysfrom1July[i],area=df$Area[i])
	   else
		   df$cf[i]=average_cf(island=island,year=year,ndays2=df$daysfrom1July[i],lev=paste(df$Position[i],df$Substrate[i],sep=""))
	}
	df$estimate.dead=df$cumdead*df$cf	
#   sum over strata for each occasion
	bystrata=df
	df=data.frame(occasion=sort(unique(df$Occasion)),cumdead=tapply(df$cumdead,df$Occasion,sum),estimate.dead=tapply(df$estimate.dead,df$Occasion,sum))
	df$cf=df$estimate.dead/df$cumdead	
	return(list(bystrata=bystrata,total=df))
}
average_cf=function(island,year,ndays1=Inf,ndays2,area=NULL,lev=NULL)
#
# for year>1997, use 1998 for pup production and early mortality prior to live count and 
# use average of 1998 and 2002 for early mortality after live count.
#
# for year <=1997, use mean of 1994 and 1995 for pup production and early pup mortality prior to live count 
# and use 1995, 1998 and 2002 analysis by area for corrections after live count (need to think about including 1998/2002)
#
{	
	if(is.nan(ndays2))return(0)
	if(island=="SMI")
	{
		if(year<1998)
		{
			if(ndays2<=34)
				cfyears=1994:1995
			else
				cfyears=c(1998,2002)
			area=as.character(area,1,3)
			area[area%in%c("EAC","WAC","ACV","PBS")]="SCV"
			area[area%in%c("NEP","NWP","PBP")]="PTS"
			if(!area%in%c("SCV","PTS","WCV","NWC"))stop(paste("Invalid area = ",area))
		}
		else
			if(ndays2<=34)
				cfyears=1998
			else
			    cfyears=c(1998,2002)
	} else
		cfyears=2006		
	cf=vector("numeric",length=length(cfyears))
	i=1
	for (y in cfyears)
	{	
		if(year < 1998 & y > 1997)
		{
			if(!eval(parse(text=paste("exists('",tolower(island),y,"a.popan.results')",sep=""))))
				eval(parse(text=paste("data(",tolower(island),y,"a.popan.results)",sep="")))
			xx=eval(parse(text=paste(tolower(island),y,"a.popan.results$cfdata$BiGross",sep="")))
			
		}else
		{
			if(!eval(parse(text=paste("exists('",tolower(island),y,".popan.results')",sep=""))))
				eval(parse(text=paste("data(",tolower(island),y,".popan.results)",sep="")))
			xx=eval(parse(text=paste(tolower(island),y,".popan.results$cfdata$BiGross",sep="")))
		}
		if(year<1998)
		{
			if(is.null(area)) stop("area must be provided if year<1998")
			select=xx[,"Area code"]==area
		} else
		{	
			if(is.null(lev)) stop("sub and pos (e.g. lev='AC') must be provided if year>1997")
			pos=substr(lev,1,1)
			sub=substr(lev,2,2)
			select=xx$Position==pos&xx$Substrate==sub
		}
		if(!is.infinite(ndays1))
		{
			if(ndays1<min(xx$daysfrom1July[select]))ndays2=min(xx$daysfrom1July[select])
			if(ndays1>max(xx$daysfrom1July[select]))ndays2=max(xx$daysfrom1July[select])
		}
		if(ndays2<min(xx$daysfrom1July[select]))ndays2=min(xx$daysfrom1July[select])
		if(ndays2>max(xx$daysfrom1July[select]))ndays2=max(xx$daysfrom1July[select])
		cf[i]=compute_cf(xx[select,],ndays1=ndays1,ndays2=ndays2)
		i=i+1
	}
	return(mean(cf))
}
# Constructs correction factor that is useful for either the pup production estimate or
# the early pup survivorship curve.  Uses arguments cfdata - correction factor data;
# ndays1 and ndays2 the range of days to use.  For pup production, default value of ndays1 
# should be used and ndays2 should be the day value of the live count. For early pup mortality (survivorship curve), 
# ndays1 is the day value of the survey at the beginning of the interval and ndays2 is the day value of the survey
# at the end of the interval. The correction factor is = (estimated dead at ndays2-estimated dead at ndays1)/
# (cummulative count at ndays2 - cumulative count at ndays1).  This is the same for pup production except
# for ndays1 the count and estimate are both 0.
compute_cf=function(BiGross,ndays1=Inf,ndays2)
{
	xx=BiGross
	if(is.infinite(ndays1))
	{
		counted1=0
		estdied1=0
	} else
	{
		xout=ndays1
		xout[xout<min(xx$daysfrom1July)]=min(xx$daysfrom1July)
		counted1=approx(x=xx$daysfrom1July,y=xx$cumdead,xout=xout)$y
		estdied1=approx(x=xx$daysfrom1July,y=xx$estimate,xout=xout)$y
	}
	xout=ndays2
	xout[xout<min(xx$daysfrom1July)]=min(xx$daysfrom1July)
	counted2=approx(x=xx$daysfrom1July,y=xx$cumdead,xout=xout)$y
	estdied2=approx(x=xx$daysfrom1July,y=xx$estimate,xout=xout)$y
	cf=(estdied2-estdied1)/(counted2-counted1)   			
	return(cf)
}

process_ch=function(ch,freq=NULL,all=FALSE)
#################################################################################
# process.ch - from a capture history creates vector of first and last times seen,
#              freq vector and indicator matrices used in log-likelihood calculations
#
#  Argument:
#        ch       - vector of character strings of 0/1
#        freq     - frequency of that ch; if null assumed to be 1; if <0 it
#                   signifies a loss on capture at last capture event
#        all      - if TRUE, computes all indicator matrices
#
#  Value: list with following elements
#        nocc        - number of capture occasions
#        freq        - absolute value of frequency for each ch
#        first       - vector of occasion numbers for first 1
#        last        - vector of occasion numbers for last 1
#        loc         - indicator of a loss on capture if set to 1
#        chmat       - capture history matrix
#        The following only returned if all==TRUE
#        FtoL        - 1's from first (1) to last (1) and 0's elsewhere (excluding
#        Fplus       - 1's from occasion after first (1) to nocc(last occasion)
#        Lplus       - 1's from occasion after last (1) to nocc
#        L           - 1's from last (1) to nocc
#        First       - 1's from occasion first (1) to nocc(last occasion)
#################################################################################
{
#  is ch comma separated? If not, separate by commas
	if(length(grep(",",ch[1]))==0)
		ch=sapply(strsplit(ch,""),paste,collapse=",")
	ch.lengths=sapply(strsplit(ch,","),length)
	nocc=ch.lengths[1]
	if(any(ch.lengths!=nocc))
		stop("\nCapture history length is not constant. \nch must be a character string with constant length or comma separated with constant number of elements \n")	
	nch=length(ch)
	if(is.null(freq))freq=rep(1,nch)
# in case multistate data are passed change all non-zero to 1
	chmat=matrix((unlist(strsplit(ch,","))),byrow=TRUE,ncol=nocc,nrow=nch)
	ch=apply(t(apply(splitCH(ch),1,function(x){ 
								x[x!="0"]=1 
								return(x)
							})),1,paste,collapse="")
#  create a matrix with 1:nocc in each row and one row for each ch
	nums=matrix(1:nocc,nrow=nch,ncol=nocc,byrow=TRUE)
#  store in a temp matrix and assign any 0 value to NA
	ymat=matrix(as.numeric(unlist(strsplit(ch,""))),byrow=TRUE,ncol=nocc,nrow=nch)
	if(suppressWarnings(all(is.numeric(as.numeric(chmat)))))chmat=ymat
	ymat[ymat==0]=NA
#  multiply nums matrix times the chmat
	nums=nums*ymat
#  use apply to get the minimum occasion and max occasion excluding NA values
	first=apply(nums,1,min,na.rm=TRUE)
	last=apply(nums,1,max,na.rm=TRUE)
	loc=rep(0,length(first))
	loc[freq<0]=1
	freq=abs(freq)
#  using the first and last values for each ch, compute the indicator matrices as defined above
	if(all)
	{
		FtoL=t(apply(cbind(first,last),1,function(x,nocc) {return(c(rep(0,x[1]),rep(1,x[2]-x[1]),rep(0,nocc-x[2])))},nocc=nocc))
		First=t(sapply(first,function(x,nocc){return(c(rep(0,x[1]-1),rep(1,nocc-x[1]+1)))},nocc=nocc))
		Fplus=t(sapply(first,function(x,nocc){return(c(rep(0,x[1]),rep(1,nocc-x[1])))},nocc=nocc))
		Lplus=t(sapply(last,function(x,nocc){return(c(rep(0,x[1]),rep(1,nocc-x[1])))},nocc=nocc))
		L=t(sapply(last,function(x,nocc){return(c(rep(0,x[1]-1),rep(1,nocc-x[1]+1)))},nocc=nocc))
#  return a list with each of the values
		return(list(nocc=nocc,freq=freq,first=first,last=last,loc=loc,chmat=chmat,FtoL=FtoL,Fplus=Fplus,Lplus=Lplus,L=L,First=First))
	}
	else
		return(list(nocc=nocc,freq=freq,first=first,last=last,loc=loc,chmat=chmat))      
}
splitCH <- function(x="ch", data=NULL, prefix="Time"){
#   Set value of ch depending on what arguments are set
	if(is.null(data)){
		ch=x
	} else
	{
		if(!x%in%names(data)) stop(paste("value for data does not contain field", x))
		ch=data[,x]
	}
#   split ch assuming non-numeric fields
	sep=""
	if(length(grep(",",ch[1]))!=0) sep=","
	chmat=do.call("rbind",strsplit(ch,sep))
#   if all fields are numeric split as numeric
	if(!any(!chmat%in%as.character(0:9)))
		chmat=t(sapply(strsplit(ch,sep),function(x)as.numeric(x)))
	if((is.character(x) & length(x)==1) & !is.null(data)){
		colnames(chmat) <- paste(prefix, c(1:ncol(chmat)),sep="")
		rownames(chmat) <- NULL
		return(cbind(data,chmat))
	}
	else{
		colnames(chmat) <- paste(prefix, c(1:ncol(chmat)),sep="")
		return(chmat)
	}
}

