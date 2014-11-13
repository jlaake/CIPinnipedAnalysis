#' Dead pup capture histories
#' 
#' Creates capture histories for tagged and untagged dead pups for POPAN analysis
#' 
#' Extracts data for a particular year, island and pup development stage (premie/fullterm) from 
#' Zc dead tag initial and Zc dead tag resight for tagged dead pups, 
#' and Zc Cu dead pup census for untagged stacked dead pups. Creates a capture history
#' for each pup and uses -1 for freq for any stacked pup.  Includes initial substrate, carcass
#' condition and beach position as covariates. Note FP (flood pond) has been used infrequently
#' and has been coverted to A (above beach crest) versus B (below beach crest).
#'   
#' It reports any mismatches (id in resights with no initial) and any duplicate initial data
#' records.
#' @export getdead_ch 
#' @param island ("SMI" or "SNI")
#' @param year four digit numeric year
#' @param development either "premie" or "fullterm"
#' @return dataframe containing ch(capture history), Carcass condition: F-fresh, D- decomposing, P- pancake, 
#' Position: A or B, and Substrate: N (non-consolidated - sandy), C (consolidated - rocky)
#' @author Jeff Laake 
getdead_ch=function(island,year,development="fullterm")
{
	development=tolower(development)
	island=toupper(island)
	initial=getCalcurData("CIPCensus","Zc dead tag initial")
	resight=getCalcurData("CIPCensus","Zc dead tag resight")
	deadstacked=getCalcurData("CIPCensus","Zc Cu dead pup census")
	deadpupareas=subset(getCalcurData("CIPCensus","DeadPupSampleAreas"),subset=YearAreaSpecies=="Zc",select=c("YearArea","Dead pup sample area","Location"))
	initial=initial[initial$Island==island & initial$Year==year & tolower(initial$Development)==development,]
	if(nrow(initial)==0) 
	{
		message(paste("\nno tagging data avaliable for island=",island, " and year = ",year,"\n",sep=""))	
		tags=FALSE
	}else
	{
		tags=TRUE
		resight=resight[resight$Island==island & resight$Year==year & tolower(resight$Development)==development,]
	}
	deadstacked$YearArea=paste(deadstacked$Year,deadstacked[,"Area code"],sep="")
	deadstacked=merge(deadstacked,deadpupareas,all.x=TRUE,by="YearArea")
	deadstacked=deadstacked[deadstacked$Location=="Mainland"&deadstacked["Dead pup sample area"]=="Y",]
	deadstacked=deadstacked[deadstacked$Species=="Zc"&deadstacked$Island==island & deadstacked$Year==year & tolower(deadstacked$Development)==development,]
	deadstacked=deadstacked[!is.na(deadstacked[,"Survey number"]),]
    if(tags)
	{
		if(any(is.na(initial[,"Survey date"]))) cat(paste(initial$ID[is.na(initial[,"Survey date"])],collapse="\n"))
		if(any(is.na(resight[,"Survey date"]))) cat(paste(resight$ID[is.na(resight[,"Survey date"])],collapse="\n"))
	}
	if(any(is.na(deadstacked[,"Survey date"]))) cat(paste(deadstacked$ID[is.na(deadstacked[,"Survey date"])],collapse="\n"))
	if(tags)
	{
		dates=rbind(initial[,c("Survey number","Survey date")],resight[,c("Survey number","Survey date")],deadstacked[,c("Survey number","Survey date")])
		surveynumbers=unique(c(deadstacked[,"Survey number"],initial[,"Survey number"],resight[,"Survey number"]))
	}
    else
	{
		dates=deadstacked[,c("Survey number","Survey date")]
		surveynumbers=unique(c(deadstacked[,"Survey number"]))
	}
	surveynumbers=surveynumbers[order(surveynumbers)]
	dates$days=(dates[,"Survey date"]-min(dates[,"Survey date"]))/(60*60*24)
	daysfrom1July=(as.Date(dates[,"Survey date"])- as.Date(paste(year,"-07-01",sep="")))
	daysfrom1July=floor(tapply(daysfrom1July,dates[,"Survey number"],mean)+.5)
	days=floor(tapply(dates$days,dates[,"Survey number"],mean)+.5)
	days=days-min(days)
	if(tags)
	{
		join=merge(initial,resight,by="SpecimenID",all.x=TRUE)
		bind=subset(initial,select=c("SpecimenID","Survey number","Carcass condition"))
		bind$Disposition="Left"
		bind=rbind(bind,subset(resight,select=c("SpecimenID","Survey number","Disposition","Carcass condition")))
		bind$Disposition=factor(bind$Disposition,levels=c("Left","Stacked"))
		bind$SpecimenID=factor(as.character(bind$SpecimenID))
		bind[,"Survey number"]=factor(as.character(bind[,"Survey number"]),levels=surveynumbers)
		ch=with(bind,table(list(SpecimenID,bind[,"Survey number"])))
		ch=apply(ch,1,paste,collapse="")
		if(length(grep("2",ch))>0) stop(paste("\n mutliple sightings at same survey \n",paste(names(ch)[grep("2",ch)],collapse="\n"),sep=""))
		freq=1-2*with(bind,table(list(SpecimenID,Disposition))) [,2]
		df=data.frame(ch=ch,freq=freq,stringsAsFactors=FALSE)
		init=subset(initial[order(initial$SpecimenID),],select=c("SpecimenID","Carcass condition","Position","Substrate","Area code"))
		init[,"Area code"]=as.character(substr(init[,"Area code"],1,3))
		init[,"Area code"][init[,"Area code"]%in%c("EAC","WAC","ACV","PBS")]="SCV"
		init[,"Area code"][init[,"Area code"]%in%c("NEP","NWP","PBP")]="PTS"
		init[,"Area code"]=factor(init[,"Area code"])
		if(any(!rownames(df)%in%init$SpecimenID)) cat("\nResights without matching ID\n",paste(rownames(df)[!rownames(df)%in%init$SpecimenID],collapse="\n"))
		if(any(!init$SpecimenID%in%rownames(df))) cat("\nInitial",paste(rownames(init)[!init$SpecimenID%in%rownames(df)],collapse="\n"))
		tt=table(init$SpecimenID)
		cat(paste(names(tt[tt>1]),collapse="\n"),"\n")
		if(nrow(init)>nrow(df))
			stop("duplicate initial records")
		else
		{
			if(nrow(init)<nrow(df))
			{
				stop("resight without initial record")
			}
			else
			{
				init$SpecimenID=NULL
				df1=cbind(df,init)
			}
		}
	}else
		df1=NULL
	if(any(is.na(deadstacked[,"Number dead"])))
	{
		cat("\n Following records have no value for Number dead and have been removed.\n")
		for(w in which(is.na(deadstacked[,"Number dead"])))
			cat("ID =",deadstacked$ID[w],"\n")
		deadstacked=deadstacked[!is.na(deadstacked[,"Number dead"]),]
	}
	deadstacked$ID=factor(as.character(deadstacked$ID))
	deadstacked[,"Survey number"]=factor(as.character(deadstacked[,"Survey number"]),levels=surveynumbers)
	ch=apply(table(deadstacked$ID,deadstacked[,"Survey number"]),1,paste,collapse="")
	df=data.frame(ID=names(ch),ch=ch)
	df$ch=as.character(df$ch)
	df=cbind(df,subset(deadstacked[order(deadstacked$ID),],select=c("Carcass condition","Position","Substrate","Number dead","Area code")))
	df$freq=-df[,"Number dead"]
	df=rbind(df1,subset(df,select=c("ch","freq","Carcass condition","Position","Substrate","Area code")))
	df[,"Area code"]=as.character(substr(df[,"Area code"],1,3))
	df[,"Area code"][df[,"Area code"]%in%c("EAC","WAC","ACV","PBS")]="SCV"
	df[,"Area code"][df[,"Area code"]%in%c("NEP","NWP","PBP")]="PTS"
	df[,"Area code"]=factor(df[,"Area code"])
	df$Position[df$Position=="FP"]="A"
	df$Substrate[df$Substrate%in%c("K","F")]="N"
	df$Position=factor(as.character(df$Position))
	df$Substrate=factor(as.character(df$Substrate))
	df[,"Carcass condition"]=factor(as.character(df[,"Carcass condition"]),levels=c("F","D","P"))
	return(list(df=df,days=days,daysfrom1July=daysfrom1July))
}
