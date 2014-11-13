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
#' one year of tagging data but it can be used with SNI as well. 
#' 
#' @export correct_dead average_cf 
#' @param island ("SMI" or "SNI")
#' @param year four digit numeric year
#' @param cfyear the year of tagging data used to create correction factor 
#' @param cfyears vector of years of tagging data used to create average correction factor 
#' @param ndays ndays is number of days past 1 July; if NULL provides a corrected value for each survey date (occasion); if non-null provides correction at that date
#' @return a list containing a dataframe with corrections by strata (eg group - area or substrate/position) and a dataframe with a total across strata. 
#' @author Jeff Laake 
#' @seealso popan.cf getdead_ch
correct_dead=function(island,year,cfyear,cfdata,ndays=NULL)
{
	x.popan=suppressMessages(getdead_ch(island,year))
	chlist=CIPinnipedAnalysis:::process_ch(x.popan$df$ch,x.popan$df$freq)
	times=x.popan$days
	time.intervals=as.numeric(diff(times))
	if(cfyear>1997)
		x.proc=process.data(x.popan$df,model="POPAN",time.intervals=time.intervals,groups=c("Position","Substrate"),begin.time=0)
	else
		x.proc=process.data(x.popan$df,model="POPAN",time.intervals=time.intervals,groups=c("Area code"),begin.time=0)
	
	numdead= tapply(abs(x.proc$data$freq),list(factor(chlist$first,1:nrow(times)),x.proc$data$group),sum)
	numdead[is.na(numdead)]=0
	
	cumdead=as.vector(apply(t(numdead), 1, cumsum))
	k=length(times)
	gindex=rep(1:nrow(x.proc$group.covariates),each=k)
	df=x.proc$group.covariates[gindex,,drop=FALSE]
	df$occasion=rep(1:k,nrow(x.proc$group.covariates))
	df$daysfrom1July=rep(x.popan$daysfrom1July,nrow(x.proc$group.covariates))
	df$cumdead=cumdead
#
#   compute and apply correction factor to observed data
#
	if(cfyear<=1997)
	{
		xx=cfdata$BiGross
		cf=matrix(0,nrow=length(levels(df[,"Area code"])),ncol=max(df$occasion)) 
		for(l in 1:length(levels(xx[,"Area code"])))
		{
			lev=levels(xx[,"Area code"])[l]
			xout=df$daysfrom1July[df[,"Area code"]==lev]
			xout[xout<min(xx$daysfrom1July[xx[,"Area code"]==lev])]=min(xx$daysfrom1July[xx[,"Area code"]==lev])
			cf[l,]=approx(x=xx$daysfrom1July[xx[,"Area code"]==lev],y=xx$cf[xx[,"Area code"]==lev],xout=xout)$y
		}
		df$cf=as.vector(t(cf))
		df$estimate.dead=df$cumdead*df$cf
		if(!is.null(ndays)) 
		{
			newdf=data.frame("Area code"=levels(xx[,"Area code"]))
			for(l in 1:length(levels(xx[,"Area code"])))
			{
				lev=levels(xx[,"Area code"])[l]
				newdf$observed.dead[l]=approx(x=df$daysfrom1July[df[,"Area code"]==lev],y=df$cumdead[df[,"Area code"]==lev],xout=ndays)$y
				newdf$estimate.dead[l]=approx(x=df$daysfrom1July[df[,"Area code"]==lev],y=df$estimate.dead[df[,"Area code"]==lev],xout=ndays)$y
			}	
			df=newdf
		}				
	}else
	{
		xx=cfdata$BiGross
		cf=matrix(0,nrow=nrow(cfdata$NGross),ncol=max(df$occasion)) 
		for(l in 1:4)
		{
			lev=c("AC","BC","AN","BN")[l]
			pos=substr(lev,1,1)
			sub=substr(lev,2,2)
			select=xx$Position==pos&xx$Substrate==sub
			xout=df$daysfrom1July[df$Position==pos&df$Substrate==sub]
			xout[xout<min(xx$daysfrom1July[select])]=min(xx$daysfrom1July[select])
			cf[l,]=approx(x=xx$daysfrom1July[select],y=xx$cf[select],xout=xout)$y
		}
		df$cf=as.vector(t(cf))
		df$estimate.dead=df$cumdead*df$cf
		if(!is.null(ndays)) 
		{
			newdf=data.frame(strata=c("AC","BC","AN","BN"))
			for(l in 1:4)
			{
				lev=c("AC","BC","AN","BN")[l]
				pos=substr(lev,1,1)
				sub=substr(lev,2,2)
				select=df$Position==pos&df$Substrate==sub
				newdf$observed.dead[l]=approx(x=df$daysfrom1July[select],y=df$cumdead[select],xout=ndays)$y
				newdf$estimate.dead[l]=approx(x=df$daysfrom1July[select],y=df$estimate.dead[select],xout=ndays)$y
			}
			df=newdf
		}				
	}
	bystrata=df
	if(is.null(ndays))
	{
		df=data.frame(occasion=sort(unique(df$occasion)),cumdead=tapply(df$cumdead,df$occasion,sum),estimate.dead=tapply(df$estimate.dead,df$occasion,sum))
		df$cf=df$estimate.dead/df$cumdead	
	}else
	{
		df=data.frame(cumdead=sum(df$observed.dead),estimate.dead=sum(df$estimate.dead))
		df$cf=df$estimate.dead/df$cumdead	
	}
	return(list(bystrata=bystrata,total=df))
}
average_cf=function(island,year,ndays=NULL,cfyears)
{
	df=NULL
	for(y in cfyears)
	{
		if(!eval(parse(text=paste("exists('",tolower(island),y,".popan.results$cfdata')",sep=""))))
			eval(parse(text=paste("data(",tolower(island),y,".popan.results)",sep="")))
		cd=correct_dead(island,year,y,eval(parse(text=paste(tolower(island),y,".popan.results$cfdata",sep=""))),ndays)$total
		df=rbind(df,data.frame(year=y,observed.dead=cd$cumdead,estimate.dead=cd$estimate.dead,cf=cd$cf))
	}
	df=df[!is.na(df$cf),]
	return(df)
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


